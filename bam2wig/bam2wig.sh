#!/bin/sh
# create a coverage wig file from a sorted, duplicate marked bam file
# Jared Evans evans.jared@mayo.edu
# 11/2014

usage()
{
echo "
NAME:
bam2wig.sh 

DESCRIPTION:
This script creates a coverage WIG file from a sorted BAM file

USAGE:
./bam2wig.sh -i input.bam -o /output/directory/path/ -b 10 -t config.txt -e Exon.Key.txt

OPTIONS:
-i	input.bam	BAM file (required)
-o	/out/dir/	Output directory (required)
-b	10		Bin size. Set to 0 to turn off binning. Default 10
-m	20		Min mapping quality. Default 20
-t	config.txt	Tool config file (required)
-e	Exon.Key.txt	Exon key file
-n			Not to merge overlapping BED regions
-h			print out this help message
"
exit 1;
}

# set defaults
bin_size=10
min_mapq=20

while getopts "i:o:b:m:t:e:nh" opt; do
	case $opt in
		i) input_bam=$OPTARG;;
		o) output_dir=$OPTARG;;
		b) bin_size=$OPTARG;;
		m) min_mapq=$OPTARG;;
		t) tool_config=$OPTARG;;
		e) exon_bed=$OPTARG;;
		n) no_merge="YES";;
		h) usage;;
		\?) echo "See available options:" >&2
		usage;;
		:) echo "See available options:" >&2
		usage;;
	esac
done

if [ -z "$input_bam" -o -z "$output_dir" -o -z "$tool_config" ]
then
	echo "Missing Required Parameters!"
	usage
fi


	set -x
	echo "Starting bam2wig"
	echo $(date)
	
	filename=$(basename $input_bam)

	# get paths to local installs of tools
	script_path=$(grep "^PATTERNCNV=" $tool_config | cut -d"=" -f2)
	script_path=${script_path}/bam2wig
	samtools=$(grep "^SAMTOOLS=" $tool_config | cut -d"=" -f2)
	bedtools=$(grep "^BEDTOOLS=" $tool_config | cut -d"=" -f2)

	# check if samtools and bedtools paths exist
	samtools_check=$(echo -n $samtools | wc -m)
	bedtools_check=$(echo -n $bedtools | wc -m)
	if [ $samtools_check == 0 -o $bedtools_check == 0 ]
	then
		echo "ERROR: Please add the paths to local installs of SAMtools and BEDtools to config.txt and try again."
		exit 1
	fi

	# newer BEDtools versions don't work unless added to PATH
	export PATH=$bedtools:$PATH

	# bams can be comma seperated if split into parts
	input_bams_array=($(echo $input_bam | tr "," " "))

	# check if bam is sorted
	for bam in "${input_bams_array[@]}"; do
		sorted=$($samtools view -H $bam | grep "^@HD" | grep "SO:coordinate" | wc -l)
		if [ $sorted == 0 ]
		then
			echo "ERROR! The BAM file ${bam} isnt coordinate sorted. Please sort the BAM and try again."
			exit 1
		fi
	done

	# create dirs
	mkdir -p $output_dir

	# get chrs
	chrs=$($samtools view -H ${input_bams_array[0]} | grep "^@SQ" | awk '{split($0,sn,"SN:"); print sn[2]}' | awk 'ORS=":"{print $1}')

	if [ $exon_bed ]
	then	
		# check if exon key or bed file
		exon_key=$(head -1 $exon_bed | awk '($1 ~ "Chr" && $4 ~ "Bin_Count"){print $0}' | wc -l)
		if [ $exon_key -gt 0 ]
		then
			sed 1d $exon_bed | awk -F"\t" '{print $1"\t"$2-1"\t"$3"\t"$4}' > $output_dir/$filename.exon.tmp.bed
			exon_key_file=$exon_bed
			exon_bed=$output_dir/$filename.exon.tmp.bed
		fi

		column_count=$(head -1 $exon_bed | awk '{print NF}')
		name_col=""
		if [ $column_count -gt 3 ]
		then
			name_col="-c 4 -o distinct -delim \",\""
		fi
		# sort, merge, and bin bed file
		# sorting in same chr order as bam (required for coverage calculation)
		for chr in $(echo $chrs | sed 's/:/ /g')
		do
			if [ $no_merge ] # option to not merge overlapping BED regions
			then
				grep -w "^${chr}" $exon_bed | sort -k2n -T $output_dir | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' >> $output_dir/$filename.exons.bed
			else
				row_count=$(grep -w "^${chr}" $exon_bed | head | wc -l)
				if [ $row_count -gt 0 ]
				then
					grep -w "^${chr}" $exon_bed | sort -k2n -T $output_dir | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' | $bedtools/mergeBed -i stdin -d -1 $name_col >> $output_dir/$filename.exons.bed
				fi
			fi
		done
		if [ -f $output_dir/$filename.exon.tmp.bed ]
		then
			rm $output_dir/$filename.exon.tmp.bed
		fi
		if [ $bin_size == 0 ] # if bin size is 0 don't bin
		then
			cp $output_dir/$filename.exons.bed $output_dir/$filename.exons_binned.bed
		else
			#$script_path/bin_exons.pl $output_dir/$filename.exons.bed $output_dir/$filename.exons_binned.bed $bin_size
			# bin BED regions. If there are leftover lengths shorter than the bin_size then they are ignored
			cat $output_dir/$filename.exons.bed | perl -slane '$gene=""; $gene="\t".$F[3] if defined $F[3]; $length=($F[2]-$F[1]); while($length >= $bin_size){print $F[0]."\t".$F[1]."\t".($F[1]+$bin_size).$gene; $F[1]+=$bin_size; $length-=$bin_size;};' -- -bin_size=$bin_size > $output_dir/$filename.exons_binned.bed
		fi
	else
		# whole genome
		echo "todo: create bins over whole genome instead of exons"
	fi

	# calculate coverage, filter duplicates and low qual reads
	echo "Calculating coverage from pileup..."
	input_bams_list="${input_bams_array[@]}"
	if [ ${#input_bams_array[@]} -gt 1 ]
	then
		$samtools cat $input_bams_list | $samtools mpileup -d 10000000 -Q 0 -q $min_mapq - | $script_path/pileup_coverage.pl $output_dir/$filename.exons_binned.bed $output_dir/$filename.coverage.txt $chrs
	else
		$samtools mpileup -d 10000000 -Q 0 -q $min_mapq $input_bams_list | $script_path/pileup_coverage.pl $output_dir/$filename.exons_binned.bed $output_dir/$filename.coverage.txt $chrs
	fi
	rm $output_dir/$filename.exons_binned.bed

	# sort bed to match exon key order (very important)
	cat $output_dir/$filename.coverage.txt | sort -k1,1 -k2,2g -T $output_dir > $output_dir/$filename.coverage.sort.txt
	cat $output_dir/$filename.exons.bed | sort -k1,1 -k2,2g -T $output_dir > $output_dir/$filename.exons.sort.bed
	rm $output_dir/$filename.exons.bed

	# convert bed to wig
	coverage_col=$(head -1 $output_dir/$filename.coverage.sort.txt | awk '{print NF}')
	$script_path/bed2wig.pl $output_dir/$filename.coverage.sort.txt $output_dir/$filename.exons.sort.bed $output_dir/$filename.coverage.wig $bin_size $coverage_col
	rm $output_dir/$filename.exons.sort.bed
	rm $output_dir/$filename.coverage.txt
	rm $output_dir/$filename.coverage.sort.txt
	
	# quick QC check on wig file
	#$script_path/count_wig_chrs.pl $output_dir/$filename.coverage.wig | awk -F"\t" '($2 == 0){print "WARNING! There is 0 coverage in the wig file for this chromosome: "$0}'
	# Give a warning if any chromosome has 0 coverage
	cat $output_dir/$filename.coverage.wig | perl -F/\\s/ -slane 'BEGIN{$current_chr=""; $current_chr_count=0;}; if($F[0] eq "fixedStep"){@chr=split("=",$F[1]); if($chr[1] ne $current_chr){print "WARNING! There is 0 coverage in the $filename.coverage.wig file for chromosome: ".$current_chr if $current_chr ne "" and $current_chr_count == 0; $current_chr=$chr[1]; $current_chr_count=0;}}else{$current_chr_count+=$F[0];};END{print "WARNING! There is 0 coverage in the $filename.coverage.wig file for chromosome: ".$current_chr if $current_chr ne "" and $current_chr_count == 0;};' -- -filename=$filename

	if [ $exon_key -gt 0 ]
	then
		wig_lines=$(cat $output_dir/$filename.coverage.wig | wc -l)
		exon_key_lines=$(awk -F"\t" '{sum+=($4+1)}END{print sum-1}' $exon_key_file)
		if [ $wig_lines != $exon_key_lines ]
		then
			echo "ERROR! wig line count (${wig_lines}) doesn't match expected count from exon key (${exon_key_lines}) in ${filename}.coverage.wig"
		fi
	fi
	gzip -f $output_dir/$filename.coverage.wig

	echo "Finished creating wig file"
	echo $(date)



#!/bin/sh
# create a coverage wig file from a sorted, duplicate marked bam file
# Jared Evans evans.jared@mayo.edu
# 2/21/12

if [ $# -lt 5 ];
then
	echo "USAGE: ./bam2wig.sh <input bam file> <output dir> <bin size (10)> <min mapping quality (20)> <tool config.txt> <exon key>";
else					
	#set -x
	echo "Starting bam2wig"
	echo $(date)
	input_bam=$1
	output_dir=$2
	bin_size=$3
	min_mapq=$4
	tool_config=$5
	exon_bed=$6
	
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

	# check if bam is sorted
	sorted=$($samtools view -H $input_bam | grep "^@HD" | grep "SO:coordinate" | wc -l)
	if [ $sorted == 0 ]
	then
		echo "ERROR: The BAM file isnt coordinate sorted. Please sort the BAM and try again."
		exit 1
	fi

	# create dirs
	mkdir -p $output_dir

	# get chrs
	chrs=$($samtools view -H $input_bam | grep "^@SQ" | awk '{split($0,sn,"SN:"); print sn[2]}' | awk 'ORS=":"{print $1}')

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
			name_col="-nms"
		fi
		# sort, merge, and bin bed file
		# sorting in same chr order as bam (required for coverage calculation)
		for chr in $(echo $chrs | sed 's/:/ /g')
		do
			grep -w "^${chr}" $exon_bed | sort -k2n | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' | $bedtools/mergeBed $name_col >> $output_dir/$filename.exons.bed
		done
		if [ -f $output_dir/$filename.exon.tmp.bed ]
		then
			rm $output_dir/$filename.exon.tmp.bed
		fi
		$script_path/bin_exons.pl $output_dir/$filename.exons.bed $output_dir/$filename.exons_binned.bed $bin_size
		rm $output_dir/$filename.exons.bed
	else
		# whole genome
		echo "todo: create bins over whole genome instead of exons"
	fi

	# calculate coverage, filter duplicates and low qual reads
	echo "Calculating coverage from pileup..."
	$samtools mpileup -Q 0 -q $min_mapq $input_bam | $script_path/pileup_coverage.pl $output_dir/$filename.exons_binned.bed $output_dir/$filename.coverage.txt $chrs
	rm $output_dir/$filename.exons_binned.bed

	# sort bed to match exon key order (very important)
	cat $output_dir/$filename.coverage.txt | sort -k1,1 -k2,2g > $output_dir/$filename.coverage.sort.txt

	# convert bed to wig
	coverage_col=$(head -1 $output_dir/$filename.coverage.sort.txt | awk '{print NF}')
	$script_path/bed2wig.pl $output_dir/$filename.coverage.sort.txt $output_dir/$filename.coverage.wig $bin_size $coverage_col
	rm $output_dir/$filename.coverage.txt
	rm $output_dir/$filename.coverage.sort.txt
	
	# quick QC check on wig file
	$script_path/count_wig_chrs.pl $output_dir/$filename.coverage.wig | awk -F"\t" '($2 == 0){print "WARNING! There is 0 coverage in the wig file for this chromosome: "$0}'
	if [ $exon_key -gt 0 ]
	then
		wig_lines=$(cat $output_dir/$filename.coverage.wig | wc -l)
		exon_key_lines=$(awk -F"\t" '{sum+=($4+1)}END{print sum-1}' $exon_key_file)
		if [ $wig_lines != $exon_key_lines ]
		then
			echo "ERROR! wig line count (${wig_lines}) doesn't match expected count from exon key (${exon_key_lines})"
		fi
	fi

	echo "Finished creating wig file"
	echo $(date)
fi	


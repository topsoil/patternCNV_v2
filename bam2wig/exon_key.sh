#!/bin/sh
# create exon key used in PatternCNV
# Jared Evans
# 8/2014

usage()
{
echo "
NAME:
exon_key.sh 

DESCRIPTION:
This script creates a coverage WIG file from a sorted BAM file

USAGE:
./bam2wig.sh -i input.bam -o /output/directory/path/ -b 10 -t config.txt -e Exon.Key.txt

OPTIONS:
-e	exon.bed	Exon BED file (required)
-c	capturekit.bed	Capture Kit BED file (required)
-o	out.txt		Output file (required)
-b	10		Bin size. Set to 0 to turn off binning. Default 10
-x	100		Amount to extend starts and stops of regions in BED file. Default 100
-s	1000		Split regions larger than this size. Default 1000
-t	config.txt	Tool config file (required)
-n			Not to merge overlapping BED regions
-h			print out this help message
"
exit 1;
}

while getopts "e:c:o:b:s:x:t:nh" opt; do
	case $opt in
		e) exon_bed=$OPTARG;;
		c) capture_bed=$OPTARG;;
		o) output_file=$OPTARG;;
		b) bin_size=$OPTARG;;
		x) extend_size=$OPTARG;;
		s) split_size=$OPTARG;;
		t) tool_config=$OPTARG;;
		n) no_merge="YES";;
		h) usage;;
		\?) echo "See available options:" >&2
		usage;;
		:) echo "See available options:" >&2
		usage;;
	esac
done

if [ -z "$exon_bed" -o -z "$capture_bed" -o -z "$output_file" -o -z "$tool_config" ]
then
	echo "Missing Required Parameters!"
	usage
fi

	set -x
	echo "Starting creation of exon key..."
	echo $(date)
	
	# set defaults
	if [ -z "$bin_size" ]
	then
		bin_size=10
	fi
	if [ -z "$extend_size" ]
	then
		extend_size=100
	fi
	if [ -z "$split_size" ]
	then
		split_size=1000
	fi

	# get paths to local installs of tools
#	script_path=$(grep "^PATTERNCNV=" $tool_config | cut -d"=" -f2)
#	script_path=${script_path}/bam2wig
	samtools=$(grep "^SAMTOOLS=" $tool_config | cut -d"=" -f2)
	bedtools=$(grep "^BEDTOOLS=" $tool_config | cut -d"=" -f2)
	genome_size=$(grep "^GENOME_SIZE=" $tool_config | cut -d"=" -f2)

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

#	column_count=$(head -1 $exon_bed | awk '{print NF}')
#	name_col=""
#	if [ $column_count -gt 3 ]
#	then
#		name_col="-nms"
#	fi
	
	echo "Sorting exon bed file..."
	# union exon+capture beds, extend, sort, merge
	if [ $no_merge ] # option to not merge overlapping BED regions
	then
		cat $exon_bed | awk -F"\t" '{if($4 != ""){print $1"\t"$2"\t"$3"\t"$4}else{print $1"\t"$2"\t"$3"\ttmp_name"}}' | $bedtools/slopBed -b $extend_size -g $genome_size | sort -k1,1 -k2,2g | awk -F"\t" -v bin=$bin_size '{if(($3-$2) < bin){print $1"\t"$2"\t"($2+bin)"\t"$4}else{print $1"\t"$2"\t"$3"\t"$4}}' > $output_file.all_exons_merged.tmp.bed
	else
		cat $exon_bed $capture_bed | awk -F"\t" '{if($4 != ""){print $1"\t"$2"\t"$3"\t"$4}else{print $1"\t"$2"\t"$3"\ttmp_name"}}' | $bedtools/slopBed -b $extend_size -g $genome_size | sort -k1,1 -k2,2g | awk -F"\t" -v bin=$bin_size '{if(($3-$2) < bin){print $1"\t"$2"\t"($2+bin)"\t"$4}else{print $1"\t"$2"\t"$3"\t"$4}}' | $bedtools/mergeBed -d -1 -nms > $output_file.all_exons_merged.tmp.bed
	fi

#	if [ $column_count -lt 4 ]
#	then
#		awk '{print $0"\t."}' $output_file.all_exons_merged.tmp.bed > $output_file.all_exons_merged.tmp.bed.tmp
#		mv $output_file.all_exons_merged.tmp.bed.tmp $output_file.all_exons_merged.tmp.bed
#	fi

	# split large regions
	cat $output_file.all_exons_merged.tmp.bed | perl -slane 'use POSIX; $split_num=ceil(($F[2]-$F[1])/$sp_size); $split_length=($F[2]-$F[1])/$split_num; for($i=0; $i<$split_num-1; $i++){print $F[0]."\t".ceil($F[1]+($split_length*$i))."\t".ceil($F[1]+($split_length*($i+1)))."\t".$F[3]}; print $F[0]."\t".ceil($F[1]+($split_length*($split_num-1)))."\t".$F[2]."\t".$F[3]' -- -sp_size=$split_size > $output_file.all_exons_merged.tmp.bed.tmp
	mv $output_file.all_exons_merged.tmp.bed.tmp $output_file.all_exons_merged.tmp.bed

	# intersect with capture kit
	$bedtools/intersectBed -a $output_file.all_exons_merged.tmp.bed -b $capture_bed -c > $output_file.all_exons_merged.incapture.tmp.bed
	rm $output_file.all_exons_merged.tmp.bed
	
	# create final exon key
#	$script_path/exon_key.pl $output_file.all_exons_merged.incapture.tmp.bed $output_file $bin_size
	cat $output_file.all_exons_merged.incapture.tmp.bed | perl -slane 'BEGIN{print "Chr\tStart\tStop\tBin_Count\tGenes\tInCapture"}; $leftover=0; $leftover=($F[2]-$F[1])%$bin_size if $bin_size != 0; $stop=$F[2]-$leftover; $bin_count=1; $bin_count=int(($stop-$F[1])/$bin_size) if $bin_size != 0; %uniq_genes=(); foreach $gene (split(";",$F[3])){$uniq_genes{$gene}=1 if $gene ne "tmp_name"}; if(!%uniq_genes){$uniq_genes{$F[0]."_".($F[1]+1)}=1}; $F[4]=1 if $F[4] > 1; print $F[0]."\t".($F[1]+1)."\t".$stop."\t".$bin_count."\t".join(":",(sort keys %uniq_genes))."\t".$F[4];' -- -bin_size=$bin_size > $output_file
	rm $output_file.all_exons_merged.incapture.tmp.bed
	
	echo "Finished creating exon key."
	echo $(date)



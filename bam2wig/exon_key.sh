#!/bin/sh
# create exon key used in PatternCNV
# Jared Evans
# 5/7/13

if [ $# != 5 ];
then
	echo "usage: ./exon_key.sh <exon bed file> <capture kit bed file> <output file> <bin size (10)> <tool config.txt>";
else					
	#set -x
	echo "Starting creation of exon key"
	echo $(date)
	exon_bed=$1
	capture_bed=$2
	output_file=$3
	bin_size=$4
	tool_config=$5
	
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

	column_count=$(head -1 $exon_bed | awk '{print NF}')
	name_col=""
	if [ $column_count -gt 3 ]
	then
		name_col="-nms"
	fi
	
	echo "Sorting exon bed file..."
	# sort, merge exons
	cat $exon_bed | sort -k1,1 -k2,2g | awk -F"\t" -v bin=$bin_size '{if(($3-$2) < bin){print $1"\t"$2"\t"($2+bin)"\t"$4}else{print $1"\t"$2"\t"$3"\t"$4}}' | $bedtools/mergeBed -d -1 $name_col > $output_file.all_exons_merged.tmp.bed

	if [ $column_count -lt 4 ]
	then
		awk '{print $0"\t."}' $output_file.all_exons_merged.tmp.bed > $output_file.all_exons_merged.tmp.bed.tmp
		mv $output_file.all_exons_merged.tmp.bed.tmp $output_file.all_exons_merged.tmp.bed
	fi

	# intersect with capture kit
	$bedtools/intersectBed -a $output_file.all_exons_merged.tmp.bed -b $capture_bed -c > $output_file.all_exons_merged.incapture.tmp.bed
	rm $output_file.all_exons_merged.tmp.bed
	
	# create exon key
	$script_path/exon_key.pl $output_file.all_exons_merged.incapture.tmp.bed $output_file $bin_size
	rm $output_file.all_exons_merged.incapture.tmp.bed
	
	echo "Finished creating exon key"
	echo $(date)
fi	


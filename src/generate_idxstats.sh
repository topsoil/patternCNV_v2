#!/bin/bash
# script to generate idxstats output for coordinate-sorted bam file
# 05/10/2018

usage()
{
echo "
NAME:
generate_idxstats.sh

DESCRIPTION:
This script indexes the supplied coordinate-sorted bam file if the index is unavailable and runs samtools' idxstats utility to retrieve and print stats in the index file 

USAGE:
./generate_idxstats.sh -c config.txt -i input.bam -o /output/directory/path/

OPTIONS:
-c	config.txt	config file (required)
-i  	input.bam	BAM file (required)
-o  	/out/dir/	Output directory (required)
-g  			Wait in grid on failure. i.e. set state to Eqw
-h			print out this help message
"
exit 1;
}


# set defaults
EXIT_CODE=1


while getopts "c:i:o:hg" opt; do
	case $opt in
		c) tool_config=$OPTARG;;
		i) bam=$OPTARG;;
		o) output_dir=$OPTARG;;
	        g) EXIT_CODE=100;;
		h) usage;;
		\?) echo "See available options:" >&2
		usage;;
		:) echo "See available options:" >&2
		usage;;
	esac
done


if [ -z "$bam" -o -z "$output_dir" -o -z "$tool_config" ]
then
	echo "Missing Required Parameters!"
	usage
fi


# get paths to local installs of tools
samtools=$(grep "^SAMTOOLS=" $tool_config | cut -d"=" -f2)


# check if samtools and bedtools paths exist
samtools_check=$(echo -n $samtools | wc -m)
if [ $samtools_check == 0 ]
then
	echo "ERROR: Please add the paths to local installs of SAMtools to config.txt and try again."
	exit 1
fi


# check if bam is sorted
sorted=$($samtools view -H $bam | grep "^@HD" | grep "SO:coordinate" | wc -l)
if [ $sorted == 0 ]
then
	echo "ERROR! The BAM file ${bam} isnt coordinate sorted. Please sort the BAM and try again."
	exit 1
fi


# index bam if index absent
remove_index=0
if [[ ! -f ${bam}.bai ]]
then
	echo "Can't find ${bam}.bai. Creating an index..."
	remove_index=1
	$samtools index $bam
fi


# dump idxstats output to <out/dir>/idxstats/<sample>.bam.idxstats
mkdir -p ${output_dir}/idxstats
bamfile=$(basename $bam)
idxstatsfile=${output_dir}/idxstats/${bamfile}.idxstats
echo -e "ref.seq.name\tref.seq.length\tnum.mapped.reads\tnum.unmapped.reads" > $idxstatsfile
$samtools idxstats $bam >> $idxstatsfile


# clean up index file if there was none to start with
if [ $remove_index == 1 ]
then
	echo "${bam}.bai wasn't originally present. Cleaning it up..."
	rm ${bam}.bai
fi


if [[ ! -f $idxstatsfile ]]
then
	echo "WARNING: ${idxstatsfile} wasn't created"
else
	echo "Finished creating idxstats file"
	echo $(date)
fi


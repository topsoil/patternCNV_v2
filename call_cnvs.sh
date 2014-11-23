#!/bin/sh
# call CNVs with PatternCNV
# 11/2014

usage()
{
echo "
NAME:
call_cnvs.sh

DESCRIPTION:
This script calls Germline or Somatic CNVs with PatternCNV

USAGE:
./call_cnvs.sh -c config.txt

OPTIONS:
-c	config.txt	config file (required)
-v			Verbose logs for debugging
-h			print out this help message
"
exit 1;
}

# default parameters
verbose="NO"

while getopts "c:vh" opt; do
	case $opt in
		c) config=$OPTARG;;
		v) verbose="YES";;
		h) usage;;
		\?) echo "See available options:" >&2
		usage;;
		:) echo "See available options:" >&2
		usage;;
	esac
done

if [ -z "$config" ]
then
	echo "Missing Required Parameters!"
	usage
fi

if [ $verbose == "YES" ]
then
	set -x
fi
echo "Start call_cnvs.sh"
echo $(date)

# parse configs
output_dir=$( cat $config | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
R=$( cat $config | grep -w '^R' | cut -d '=' -f2)
sample_info=$( cat $config | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)

patterncnv_path=$( cat $config | grep -w '^PATTERNCNV' | cut -d '=' -f2)

# check to see if Somatic calling is needed, otherwise just perform Germline calling
sample_type=$(awk -F"\t" '($3 == "Somatic"){print "Somatic"}' $sample_info )

if [ "$sample_type" ]
then
	echo "Performing Somatic CNV calling"
	$R/Rscript --vanilla $patterncnv_path/somatic.cnv.R $output_dir $patterncnv_path/Rlib $output_dir/configs/config.ini
else
	echo "Performing Germline CNV calling"
	$R/Rscript --vanilla $patterncnv_path/germline.cnv.R $output_dir $patterncnv_path/Rlib $output_dir/configs/config.ini
fi


echo "End call_cnvs.sh"
echo $(date)


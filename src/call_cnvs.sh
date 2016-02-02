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

# set new TMPDIR variable so R can start even if /tmp/ is full
export TMPDIR=$output_dir

# check to see if Somatic calling is needed, otherwise just perform Germline calling
sample_type=$(awk -F"\t" '($3 == "Somatic"){print "Somatic"}' $sample_info )

if [ "$sample_type" ]
then
	echo "Performing Somatic CNV calling"
	$R/Rscript --vanilla $patterncnv_path/src/somatic.cnv.R $output_dir $patterncnv_path/Rlib $output_dir/configs/config.ini
else
	echo "Performing Germline CNV calling"
	$R/Rscript --vanilla $patterncnv_path/src/germline.cnv.R $output_dir $patterncnv_path/Rlib $output_dir/configs/config.ini
fi

# QC check, if files are missing then send email and exit 100
cnv_txt_count=$(ls $output_dir/cnv-txt/*.txt | grep -v CNV_matrix | grep -v pval_matrix | wc -l)
cnv_plot_count=$(ls $output_dir/cnv-plot/*.png| wc -l)
if [ $cnv_txt_count -le 0 -o $cnv_plot_count -le 0 ]
then
	mailx=$(which mailx)
	if [ "$mailx" != "" ] ; then
		email=$(finger $USER | awk -F ';' '{print $2}' | head -n1)
		TMPDIR=$output_dir
		SUB="PatternCNV Error in CALL_CNVS script"
		MES="PatternCNV ERROR! After CNV calling the number of CNV txt files (${cnv_txt_count}) or CNV plot files (${cnv_plot_count}) are not correct.\n\nSGE Log files:\n$SGE_STDERR_PATH\n$SGE_STDOUT_PATH\n\nFiles to check:\nOutput Directory:\n$output_dir\n"
		echo -e "$MES" | mailx -s "$SUB" "$email"
		sleep 15s
	fi
	exit 100;
fi


echo "End call_cnvs.sh"
echo $(date)


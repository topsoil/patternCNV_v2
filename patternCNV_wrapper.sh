#!/bin/sh
# wraper script to run PatternCNV
# 11/2014

usage()
{
echo "
NAME:
patternCNV_wrapper.sh

DESCRIPTION:
This wrapper script runs PatternCNV

USAGE:
./patternCNV_wrapper.sh -c config.txt > patternCNV.log 2>&1

OPTIONS:
-c      config.txt      config file (required)
-b      10              Bin size used to determine resolution of coverage. Default: 10
-m      20              Minimum mapping quality of reads to use for coverage calculation. Default: 20
-z      1000            Size of exons that should be split into multiple exons. Default: 1000
-x      100             Extension buffer size defining size each exon should be extended on both sides. Default: 100
-n                      Don't merge overlapping exons. Default: Merge
-v                      Verbose logs for debugging
-j                      Optional comma separated list of job IDs for PatternCNV to wait on (qsub -hold_jid).
-w                      Optional prefix to use in job names (qsub -N)
-u                      Optional suffix to use in job names (qsub -N)
-h                      print out this help message
"
exit 1;
}

# default parameters if not changed by user
bin_size=10
min_mapping_qual=20
merge_overlaps="YES"
split_size=1000
extension_buffer=100

while getopts "c:bmxznvj:w:u:h" opt; do
	case $opt in
		c) config=$OPTARG;;
		b) bin_size=$OPTARG;;
		m) min_mapping_qual=$OPTARG;;
		x) extension_buffer=$OPTARG;;
		n) merge_overlaps="NO";;
		z) split_size=$OPTARG;;
		v) verbose="YES";;
		j) jobs_to_hold_for=$OPTARG;;
		w) job_name_prefix=$OPTARG;;
		u) job_name_suffix=$OPTARG;;
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

if [ "$verbose" ]
then
	set -x
fi
echo "Start PatternCNV Wrapper"
echo $(date)

# parse configs
output_dir=$( cat $config | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
sample_info=$( cat $config | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
capture_bed=$( cat $config | grep -w '^CAPTUREKIT_BED' | cut -d '=' -f2)
exon_bed=$( cat $config | grep -w '^EXON_BED' | cut -d '=' -f2)

#bin_size=$( cat $config | grep -w '^BIN_SIZE' | cut -d '=' -f2)
#min_mapping_qual=$( cat $config | grep -w '^MIN_MAPPING_QUAL' | cut -d '=' -f2)
#merge_overlaps=$( cat $config | grep -w '^MERGE_OVERLAPPING_REGIONS' | cut -d '=' -f2)
#split_size=$( cat $config | grep -w '^SPLIT_EXON_SIZE' | cut -d '=' -f2)
#extension_buffer=$( cat $config | grep -w '^EXTENSION_BUFFER' | cut -d '=' -f2)

patterncnv_path=$( cat $config | grep -w '^PATTERNCNV' | cut -d '=' -f2)
samtools_path=$( cat $config | grep -w '^SAMTOOLS' | cut -d '=' -f2)
bedtools_path=$( cat $config | grep -w '^BEDTOOLS' | cut -d '=' -f2)

email=$( cat $config | grep -w '^EMAIL' | cut -d '=' -f2)
queue=$( cat $config | grep -w '^QUEUE' | cut -d '=' -f2)
memory_exonkey=$( cat $config | grep -w '^QSUB_EXONKEY_MEMORY' | awk -F 'QSUB_EXONKEY_MEMORY=' '{print $2}')
memory_bam2wig=$( cat $config | grep -w '^QSUB_BAM2WIG_MEMORY' | awk -F 'QSUB_BAM2WIG_MEMORY=' '{print $2}')
memory_callcnvs=$( cat $config | grep -w '^QSUB_CALLCNVS_MEMORY' | awk -F 'QSUB_CALLCNVS_MEMORY=' '{print $2}')

# create output dirs
mkdir -p $output_dir
mkdir -p $output_dir/logs
mkdir -p $output_dir/wigs
mkdir -p $output_dir/configs
#mkdir -p $output_dir/qc
mkdir -p $output_dir/cnv-txt
mkdir -p $output_dir/cnv-plot

# copy configs to ouput dir
cat $config | grep -v "^EXON_KEY" | awk -F"=" -v out=$output_dir '{if($1 == "SAMPLE_INFO"){print "SAMPLE_INFO="out"/configs/sample_info.txt\nEXON_KEY="out"/configs/exon_key.txt"}else{print}}' > $output_dir/configs/config.txt
# add wigs to sample_info.txt
head -1 $sample_info | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tfile.name"}' > $output_dir/configs/sample_info.txt
sed 1d $sample_info | awk -v out=$output_dir '{n=split($5,a,"/"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"out"/wigs/"a[n]".coverage.wig.gz"}' >> $output_dir/configs/sample_info.txt
config=$output_dir/configs/config.txt
sample_info=$output_dir/configs/sample_info.txt
exon_key=$output_dir/configs/exon_key.txt
# create config.ini
echo "plot_output_DIR = '$output_dir/cnv-plot/'" > $output_dir/configs/config.ini
echo "txt_output_DIR = '$output_dir/cnv-txt/'" >> $output_dir/configs/config.ini
echo "exon_key_file = '$exon_key'" >> $output_dir/configs/config.ini
echo "sample_info_file = '$sample_info'" >> $output_dir/configs/config.ini

merge_param=""
if [ $merge_overlaps == "NO" ]
then
	merge_param="-n"
fi

previous_jobids=""
if [ "$jobs_to_hold_for" ]
then
	previous_jobids="-hold_jid $jobs_to_hold_for"
fi

job_name="PATTERNCNV"
if [ "$job_name_prefix" ]
then
	job_name="${job_name_prefix}.PATTERNCNV"
fi

job_suffix=""
if [ "$job_name_suffix" ]
then
	job_suffix=".${job_name_suffix}"
fi

# create exon key
EXONKEY=$(qsub -wd $output_dir/logs -q $queue -m a -M $email $memory_exonkey $previous_jobids -N $job_name.exon_key${job_suffix} $patterncnv_path/bam2wig/exon_key.sh -e $exon_bed -c $capture_bed -o $exon_key -b $bin_size -t $config -x $extension_buffer -s $split_size $merge_param)
jobid_exonkey=$(echo $EXONKEY | cut -d ' ' -f3)

# bam2wig
jobid_bam2wig=""
for bam in $(cut -f5 $sample_info | sed 1d)
do
	BAM2WIG=$(qsub -wd $output_dir/logs -q $queue -m a -M $email $memory_bam2wig -hold_jid $jobid_exonkey -N $job_name.bam2wig${job_suffix} $patterncnv_path/bam2wig/bam2wig.sh -i $bam -o $output_dir/wigs -b $bin_size -m $min_mapping_qual -t $config -e $exon_key $merge_param)
	jobid=$(echo $BAM2WIG | cut -d ' ' -f3)
	jobid_bam2wig="${jobid},${jobid_bam2wig}"
done


# call CNVs
qsub -wd $output_dir/logs -q $queue -m a -M $email $memory_callcnvs -hold_jid $jobid_bam2wig -N $job_name.call_cnvs${job_suffix} $patterncnv_path/call_cnvs.sh -c $config -v



echo "End PatterCNV Wrapper"
echo $(date)


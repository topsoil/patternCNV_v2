#!/bin/sh
# wrapper script to run PatternCNV
# 12/2014

usage()
{
echo "
NAME:
patternCNV_wrapper.sh

DESCRIPTION:
PatternCNV is designed to detect Germline or Somatic CNVs in exome or custom capture samples. This wrapper script will run 
all of the PatternCNV components, including bam2wig, QC, GC correction, and CNV calling.

USAGE:
./patternCNV_wrapper.sh -c config.txt > patternCNV.log 2>&1

OPTIONS:
-c      config.txt      config file (required)
-b      10              Bin size used to determine resolution of coverage. Default: 10
-m      20              Minimum mapping quality of reads to use for coverage calculation. Default: 20
-z      1000            Size of exons that should be split into multiple exons. Default: 1000
-x      100             Extension buffer size defining size each exon should be extended on both sides. Default: 100
-s                      Run patternCNV in serial mode, not submitting any jobs to an sge cluster. Default: Parallel mode
-n                      Don't merge overlapping exons. Default: Merge
-v                      Verbose logs for debugging
-j                      Optional comma separated list of job IDs for PatternCNV to wait on (qsub -hold_jid).
-w                      Optional prefix to use in job names (qsub -N)
-u                      Optional suffix to use in job names (qsub -N)
-l                      Optional path to logs folder if user doesn't want to use the default location
-i                      Incremental mode, does not generate wigs or exon_key if they exist.
-d			Debug mode. Keeps exon bed and bam2wig intermediate files.
-r                      Use all samples to determine baseline for chrX. Default: Use only females to determine chrX baseline if available; use males if females unavailable
-g  Wait in grid on failure. i.e. set state to Eqw
-h                      print out this help message
"
exit 1;
}

# default parameters if not changed by user
bin_size=10
min_mapping_qual=20
serial="NO"
merge_overlaps="YES"
split_size=1000
extension_buffer=100
incremental="NO"
rescaleX_useAllSamples="NO"
EXIT_CODE_OPT=""

while getopts "c:b:m:x:z:snvdij:w:u:l:r:gh" opt; do
	case $opt in
		c) config=$OPTARG;;
		b) bin_size=$OPTARG;;
		m) min_mapping_qual=$OPTARG;;
		x) extension_buffer=$OPTARG;;
                s) serial="YES";;
		n) merge_overlaps="NO";;
		z) split_size=$OPTARG;;
		v) verbose="YES";;
		d) debug="YES";;
                i) incremental="YES";;
		j) jobs_to_hold_for=$OPTARG;;
		w) job_name_prefix=$OPTARG;;
		u) job_name_suffix=$OPTARG;;
		l) user_logs_output=$OPTARG;;
                r) rescaleX_useAllSamples="YES";;
        g) EXIT_CODE_OPT=" -g";;
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

patterncnv_path=$( cat $config | grep -w '^PATTERNCNV' | cut -d '=' -f2)
samtools_path=$( cat $config | grep -w '^SAMTOOLS' | cut -d '=' -f2)
bedtools_path=$( cat $config | grep -w '^BEDTOOLS' | cut -d '=' -f2)

email=$( cat $config | grep -w '^EMAIL' | cut -d '=' -f2)
queue=$( cat $config | grep -w '^QUEUE' | cut -d '=' -f2)
memory_exonkey=$( cat $config | grep -w '^QSUB_EXONKEY_MEMORY' | awk -F 'QSUB_EXONKEY_MEMORY=' '{print $2}')
memory_bam2wig=$( cat $config | grep -w '^QSUB_BAM2WIG_MEMORY' | awk -F 'QSUB_BAM2WIG_MEMORY=' '{print $2}')
memory_callcnvs=$( cat $config | grep -w '^QSUB_CALLCNVS_MEMORY' | awk -F 'QSUB_CALLCNVS_MEMORY=' '{print $2}')


QSUB=""
if [[ "${serial}" == "NO" ]] ; then
    QSUB=$( cat $config | grep -w '^QSUB' | awk -F 'QSUB=' '{print $2}')
    if [[ ${#QSUB} -eq 0 ]] ; then
        QSUB="qsub"
        echo "WARNING: In PatternCNV Wrapper, QSUB not defined in config file, trying to use one in path."
    else
        QSUBTEST=`which "${QSUB}"`
        if [[ ${#QSUBTEST} -eq 0 ]] ; then
            echo "ERROR: QSUB command ($QSUB) defined in config file, not found. Trying qsub in path"
            QSUB="qsub"
        fi
    fi
    QSUBTEST=`which "${QSUB}"`
    if [[ ${#QSUBTEST} -eq 0 ]] ; then
        echo "ERROR: PatternCNV Wrapper: Cannot find QSUB command (${QSUB}) in path, perhaps you need to run in serial mode (-s) for no cluster mode";
        echo $(date)
        echo "End PatternCNV Wrapper"
        exit 127
    fi
fi

# create output dirs
logs_dir=$output_dir/logs
if [ "$user_logs_output" ]
then
	logs_dir=$user_logs_output
fi
mkdir -p $output_dir
mkdir -p $logs_dir
mkdir -p $output_dir/wigs
mkdir -p $output_dir/configs
mkdir -p $output_dir/cnv-txt
mkdir -p $output_dir/cnv-plot


# copy configs to output dir
cat $config | grep -v "^EXON_KEY" | awk -F"=" -v out=$output_dir '{if($1 == "SAMPLE_INFO"){print "SAMPLE_INFO="out"/configs/sample_info.txt\nEXON_KEY="out"/configs/exon_key.txt"}else{print}}' > $output_dir/configs/config.txt
# add wigs to sample_info.txt
# Hugues 7/08/16 .. allow arbitrary number of fields in sample info to support sex checks.
#head -1 $sample_info | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tfile.name"}' > $output_dir/configs/sample_info.txt
head -1 $sample_info | awk '{b=$1; for(i=2;i<=NF;i+=1){b=b "\t" $i};print b "\tfile.name"}' > $output_dir/configs/sample_info.txt
#sed 1d $sample_info | awk -v out=$output_dir '{n=split($5,a,"/"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"out"/wigs/"a[n]".coverage.wig.gz"}' >> $output_dir/configs/sample_info.txt
sed 1d $sample_info | awk -v out=$output_dir '{n=split($5,a,"/");b=$1; for(i=2;i<=NF;i+=1){b=b "\t" $i};print b "\t" out "/wigs/" a[n] ".coverage.wig.gz"}' >> $output_dir/configs/sample_info.txt
config=$output_dir/configs/config.txt
sample_info=$output_dir/configs/sample_info.txt
exon_key=$output_dir/configs/exon_key.txt
# create config.ini
echo "plot_output_DIR = '$output_dir/cnv-plot/'" > $output_dir/configs/config.ini
echo "txt_output_DIR = '$output_dir/cnv-txt/'" >> $output_dir/configs/config.ini
echo "exon_key_file = '$exon_key'" >> $output_dir/configs/config.ini
echo "sample_info_file = '$sample_info'" >> $output_dir/configs/config.ini
if [ $rescaleX_useAllSamples == "YES" ]
then
        echo "rescaleX_useAllSamples = 'TRUE'" >> $output_dir/configs/config.ini
else
        echo "rescaleX_useAllSamples = 'FALSE'" >> $output_dir/configs/config.ini
fi

additional_params=""
if [ $merge_overlaps == "NO" ]
then
	additional_params="${additional_params}-n "
fi
if [ "$debug" ]
then
	additional_params="${additional_params}-d "
fi
additional_params="${additional_params}${EXIT_CODE_OPT}"

previous_jobids=""
if [[ ${#jobs_to_hold_for} -gt 0 ]] ; then
	previous_jobids="-hold_jid $jobs_to_hold_for"
fi

job_name="PATTERNCNV"
if [[ ${#job_name_prefix} -gt 0 ]] ; then
	job_name="${job_name_prefix}.PATTERNCNV"
fi

job_suffix=""
if [[ ${#job_name_suffix} -gt 0 ]] ; then
	job_suffix=".${job_name_suffix}"
fi

# create exon key
echo "Checking if need to create exon key file $exon_key"
pcnv_command="$patterncnv_path/src/bam2wig/exon_key.sh -e $exon_bed -c $capture_bed -o $exon_key -b $bin_size -t $config -x $extension_buffer -s $split_size $additional_params"
hold_exonkey=""
if [[ "$incremental" = "NO" || ( "$incremental" = "YES" && ! -f "$exon_key" ) ]] ; then 
    if [ "$serial" == "YES" ]
    then
        $pcnv_command
        echo -e "\n# PatternCNV ExonKey Job for all samples\n${pcnv_command}\n"
    else 
	qsub_command="${QSUB} -wd $logs_dir -q $queue -m a -M $email $memory_exonkey $previous_jobids -N $job_name.exon_key.allsamples${job_suffix} $pcnv_command"
	EXONKEY=$($qsub_command)
	echo -e "\n# PatternCNV ExonKey Job Submission for all samples\n${qsub_command}"
	echo -e "${EXONKEY}\n"
	jobid_exonkey=$(echo $EXONKEY | cut -d ' ' -f3)
	hold_exonkey="-hold_jid $jobid_exonkey"
    fi
fi

# bam2wig for each unique sample
jobid_bam2wig=""

for sample in $(cut -f1 $sample_info | sed 1d | sort | uniq)
do
	bam=$(grep -P "^${sample}\t" $sample_info | head -1 | cut -f5)
	bamfile=$(basename $bam)
	wigfile="$output_dir/wigs/${bamfile}.coverage.wig.gz"
	echo "Testing if need to generate $wigfile"
	if  [[ "$incremental" = "NO" || ( "$incremental" = "YES" && ! -f "$wigfile" ) ]] ; then
	    pcnv_command="$patterncnv_path/src/bam2wig/bam2wig.sh -i $bam -o $output_dir/wigs -b $bin_size -m $min_mapping_qual -t $config -e $exon_key $additional_params"
	    if [ "$serial" == "YES" ]
            then
                echo -e "# PatternCNV BAM2WIG Job for sample ${sample}\n${pcnv_command}\n"
                $pcnv_command
	    else
		qsub_command="${QSUB} -wd $logs_dir -q $queue -m a -M $email $memory_bam2wig $hold_exonkey -N $job_name.bam2wig.${sample}${job_suffix} $pcnv_command"
		BAM2WIG=$($qsub_command)
		echo -e "# PatternCNV BAM2WIG Job Submission for sample ${sample}\n${qsub_command}"
		echo -e "${BAM2WIG}\n"

		jobid=$(echo $BAM2WIG | cut -d ' ' -f3)
		jobid_bam2wig="${jobid},${jobid_bam2wig}"
	    fi
	fi
done
echo "jobid_bam2wig: $jobid_bam2wig"


# idxstats for each unique sample
jobid_idxstats=""
for sample in $(awk '{ if ($3=="Somatic") {print $1} }' $sample_info | sort | uniq)
do
	bam=$(grep -P "^${sample}\t" $sample_info | head -1 | cut -f5)
	idxstats_command="$patterncnv_path/src/generate_idxstats.sh -c $config -i $bam -o $output_dir"
	if [[ ${#jobid_bam2wig} -gt 0 ]] ; then
		qsub_command="${QSUB} -wd $logs_dir -q $queue -m a -M $email $memory_bam2wig -hold_jid $jobid_bam2wig   $hold_exonkey -N $job_name.idxstats.${sample}${job_suffix} $idxstats_command"
	else
		qsub_command="${QSUB} -wd $logs_dir -q $queue -m a -M $email $memory_bam2wig                            $hold_exonkey -N $job_name.idxstats.${sample}${job_suffix} $idxstats_command"
	fi
	IDXSTATS=$($qsub_command)
	echo -e "# PatternCNV IDXSTATS Job Submission for sample ${sample}\n${qsub_command}"
	echo -e "${IDXSTATS}\n"
	jobid=$(echo $IDXSTATS | cut -d ' ' -f3)
	jobid_idxstats="${jobid},${jobid_idxstats}"
done
echo "jobid_idxstats: $jobid_idxstats"


# call CNVs
pcnv_command="${patterncnv_path}/src/call_cnvs.sh -c $config -v ${EXIT_CODE_OPT}"
if [ "$serial" == "YES" ]
then
    $pcnv_command
    echo -e "# PatternCNV CallCNVs Job for all samples\n${pcnv_command}\n"
else
    if [[ ${#jobid_idxstats} -gt 0 ]] ; then
	qsub_command="${QSUB} -wd $logs_dir -q $queue -m a -M $email $memory_callcnvs -hold_jid $jobid_idxstats   $hold_exonkey -N $job_name.call_cnvs.allsamples${job_suffix} $pcnv_command"
    else 
	qsub_command="${QSUB} -wd $logs_dir -q $queue -m a -M $email $memory_callcnvs                             $hold_exonkey -N $job_name.call_cnvs.allsamples${job_suffix} $pcnv_command"
    fi
    echo -e "# PatternCNV CallCNVs Job Submission for all samples\n${qsub_command}"
    CALLCNVS=$($qsub_command)
    echo -e "${CALLCNVS}\n"
fi

echo "End PatternCNV Wrapper"
echo $(date)


#!/bin/sh
# SGE wrapper around bam2wig and exon key creation for a batch of bams
# 5/13/13

if [ $# != 1 ];
then
	echo "usage: ./sge_wrapper.sh <sge_config file>";
else
	set -x
	echo $(date)
	config=$1

	# get parameters
	email=$( cat $config | grep -w '^EMAIL' | cut -d '=' -f2)
	queue=$( cat $config | grep -w '^QUEUE' | cut -d '=' -f2)
	output_dir=$( cat $config | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
	bams_config=$( cat $config | grep -w '^BAMS_CONFIG' | cut -d '=' -f2)
	capture_bed=$( cat $config | grep -w '^CAPTUREKIT_BED' | cut -d '=' -f2)
	bin_size=$( cat $config | grep -w '^BIN_SIZE' | cut -d '=' -f2)
	min_mapping_qual=$( cat $config | grep -w '^MIN_MAPPING_QUAL' | cut -d '=' -f2)
	exon_bed=$( cat $config | grep -w '^EXON_BED' | cut -d '=' -f2)
	script_path=$( cat $config | grep -w '^PATTERNCNV' | cut -d '=' -f2)
	samtools_path=$( cat $config | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	bedtools_path=$( cat $config | grep -w '^BEDTOOLS' | cut -d '=' -f2)
	memory=$( cat $config | grep -w '^QSUB_MEMORY' | awk -F 'QSUB_MEMORY=' '{print $2}')
	merge_overlaps=$( cat $config | grep -w '^MERGE_OVERLAPPING_REGIONS' | cut -d '=' -f2)

	merge_param=""
	if [ $merge_overlaps == "NO" ]
	then
		merge_param="-n"
	fi

	# create dirs
	mkdir -p $output_dir/logs

	cp $config $output_dir/sge_config.txt
	config=$output_dir/sge_config.txt

	# create exon key
	EXONKEY=$(qsub -V -wd $output_dir/logs -q $queue -m a -M $email $memory $script_path/bam2wig/exon_key.sh -e $exon_bed -c $capture_bed -o $output_dir/PatternCNV.Exon.Key.txt -b $bin_size -t $config $merge_param)
	jobid_exonkey=$(echo $EXONKEY | cut -d ' ' -f3)

	# bam2wig
	for bam in $(cat $bams_config)
	do
		qsub -V -wd $output_dir/logs -q $queue -m a -M $email $memory -hold_jid $jobid_exonkey $script_path/bam2wig/bam2wig.sh -i $bam -o $output_dir -b $bin_size -m $min_mapping_qual -t $config -e $output_dir/PatternCNV.Exon.Key.txt $merge_param
	done
	echo $(date)
fi

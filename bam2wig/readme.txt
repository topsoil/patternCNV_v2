# Create the necessary WIG files to run PatternCNV using the SGE wrapper script
# Jared Evans
# 8/14

## Step 1. Create sge_config.txt and sge_bams_config.txt files:

# Example sge_config.txt file:
EMAIL=evans.jared@mayo.edu
QUEUE=1-day
OUTPUT_DIR=/data2/bsi/staff_analysis/m088341/tools/PatternCNV/sge_test

BAMS_CONFIG=/projects/bsi/bioinf_int/s112423.PatternCNV/trunk/bam2wig/sge_bams_config.txt
CAPTUREKIT_BED=/data5/bsi/refdata/exomeCapture/SureSelect_All_Exon_V4+UTRs_with_annotation_hg19_.bed
EXON_BED=/data2/bsi/staff_analysis/m088341/tools/PatternCNV/reference/exon_beds/refflat_exons_hg19.bed
GENOME_SIZE=/projects/bsi/bioinf_int/s112423.PatternCNV/trunk/data/human_37.1.genome

BIN_SIZE=10
MIN_MAPPING_QUAL=20
MERGE_OVERLAPPING_REGIONS=YES
SPLIT_EXON_SIZE=1000
EXTENSION_BUFFER=100
PATTERNCNV=/projects/bsi/bioinf_int/s112423.PatternCNV/trunk
SAMTOOLS=/projects/bsi/bictools/apps/alignment/samtools/samtools-0.1.18/samtools
BEDTOOLS=/projects/bsi/bictools/apps/misc/BEDTools/2.15.0/bin
QSUB_MEMORY=-l h_vmem=2G -l h_stack=10M


# Example sge_bams_config.txt file:
/full/path/to/bam1.bam
/full/path/to/bam2.bam
/full/path/to/bam3.bam


## Step 2. run the sge_wrapper.sh script with will submit jobs to the cluster
#to create the exon key plus WIG files for each bam:
./sge_wrapper.sh /path/to/your/sge_config.txt


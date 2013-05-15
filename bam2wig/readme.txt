# Create the necessary files to run PatternCNV
# Jared Evans
# 5/7/13



## STEP 1. Create an Exon Key:

# NOTE: start and stop of output key file are 1-base so it matches the 1-base WIG files
./exon_key.sh <exon bed file> <capture kit bed file> <output file> <bin size (10)> <tool config.txt>

# EXAMPLE USAGE:
./exon_key.sh reference/refflat_exons_hg19.bed SureSelect_All_Exon_V4+UTRs_with_annotation_hg19_.bed PatternCNV.Exon.Key.txt 10 config.txt


## STEP 2. Generate WIG files for each BAM:
./bam2wig.sh <input bam file> <output dir> <bin size (10)> <min mapping quality (20)> <tool config.txt> <exon key>

# EXAMPLE USAGE:
./bam2wig.sh LS1946_blood.sorted.bam output_dir/ 10 20 config.txt PatternCNV.Exon.Key.txt


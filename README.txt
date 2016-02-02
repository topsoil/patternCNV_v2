INTRODUCTION
------------

PatternCNV is designed to detect Germline or Somatic CNVs in exome or custom capture samples.


WEBSITE
-------

http://bioinformaticstools.mayo.edu/research/patterncnv/


PUBLICATION
-----------

Wang, C., J. M. Evans, et al. (2014). PatternCNV: a versatile tool for detecting copy number changes 
from exome sequencing data. Bioinformatics 30(18): 2678-2680.


CONTACTS
--------

Chen Wang (wang.chen@mayo.edu)
Jared Evans (evans.jared@mayo.edu)


USAGE
-----

./patternCNV_wrapper.sh [options] -c config.txt > patternCNV.log 2>&1


OPTIONS
-------

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
-d			Debug mode. Keeps exon bed and bam2wig intermediate files.
-h                      print out this help message


INSTALLATION
------------

User must install these tools on their system:
 - Samtools (http://www.htslib.org)
 - bedtools (http://bedtools.readthedocs.org)
 - R (https://www.r-project.org)
 - Perl (https://www.perl.org)

Install required R packages with the following script:
  cd patterncnv_code_dir/
  Rscript install_R_packages.R


CONFIG
------

OUTPUT_DIR=/path/to/output/dir/
SAMPLE_INFO=/path/to/sample_info.txt
CAPTUREKIT_BED=/path/to/capture_kit_regions.bed
EXON_BED=/path/to/exon_regions.bed
GENOME_SIZE=/path/to/contig/sizes/human_37.1.genome
GENOME_REF=/path/to/reference/sequence/hg19.fa

PATTERNCNV=/path/to/patterncnv/code/
SAMTOOLS=/path/to/samtools
BEDTOOLS=/path/to/bedtools/bin
R=/path/to/R/bin

EMAIL=email@email.com
QUEUE=1-day
QSUB_EXONKEY_MEMORY=-l h_vmem=4G -l h_stack=10M
QSUB_BAM2WIG_MEMORY=-l h_vmem=2G -l h_stack=10M
QSUB_CALLCNVS_MEMORY=-l h_vmem=4G -l h_stack=10M


SAMPLE_INFO_FILE
----------------

sample.name	subject.ID	sample.type	batch.ID	BAM.file
sampleid1	sampleid1	Germline	1	/path/to/sampleid1.bam
sampleid2	sampleid2	Germline	1	/path/to/sampleid2.bam
sampleid3	sampleid3	Germline	1	/path/to/sampleid3.bam

sample1_blood	pair1	Germline	1	/path/to/sample1_blood.bam
sample1_tumor	pair1	Somatic	1	/path/to/sample1_tumor.bam


CAPTUREKIT_BED
--------------

3 column capture region:

chr1	10000	10200
chr1	11000	11200
...


EXON_BED
--------

4 column exon definitions with 4th column being gene name:

chr1	24737	24891	WASH7P
chr1	29320	29370	WASH7P
chr1	34610	35174	FAM138A
chr1	34610	35174	FAM138F
...





#!/usr/bin/perl

## Calculate the average coverage of BED intervals from a piped-in pileup file
## Jared Evans 2/25/13


use strict;
use warnings;


if (scalar(@ARGV) != 3)	
{
	die ( "USAGE: pileup_coverage.pl [BED file] [output file] [ordered list of chromosomes (colon seperated)]\n" );
}


# input must be sorted in order for this algorithm to work
open BED, "$ARGV[0]" or die "opening file: $ARGV[0]";
open OUT, ">", "$ARGV[1]" or die "opening file: $ARGV[1]";
my @chrs = split(":",$ARGV[2]);

# create chr hash with ordered values. It is required to know the order of the chrs to properly walk through both files
my %chr_order = ();
for(my $i = 0; $i < scalar(@chrs); $i++){
	$chr_order{$chrs[$i]} = $i;
}

# prime the pump
my $bed_line = <BED>;
my $pileup_line = <STDIN>;
my $total_reads = 0;

while(defined $bed_line || defined $pileup_line){
	chomp $bed_line if defined $bed_line;
	chomp $pileup_line if defined $pileup_line;
	
	if(!defined $bed_line){
		# no need to go further if reached end of BED
		last
	}elsif(!defined $pileup_line){
		# reached the end of pileup, make sure to keep calculating 0 coverage for rest of BED
		print OUT "$bed_line\t";
		if($total_reads == 0){
			print OUT "$total_reads\n";
		}else{
			my @bed_array = split("\t",$bed_line);
			my $mean_cov = int($total_reads/($bed_array[2]-$bed_array[1]));
			print OUT "$mean_cov\n";
		}
		$total_reads = 0;
		$bed_line = <BED>;
	}else{
		my @bed_array = split("\t",$bed_line);
		my @pileup_array = split("\t",$pileup_line);
		if(($bed_array[0] eq $pileup_array[0]) && ($pileup_array[1] > $bed_array[1]) && ($pileup_array[1] <= $bed_array[2])){
			# pileup region falls within bed region
			$total_reads+=$pileup_array[3];
			$pileup_line = <STDIN>;
			
		}elsif((($bed_array[0] eq $pileup_array[0]) && ($pileup_array[1] > $bed_array[2])) || ($chr_order{$pileup_array[0]} > $chr_order{$bed_array[0]})){
			# pileup is now past BED interval, go ahead and print coverage
			print OUT "$bed_line\t";
			if($total_reads == 0){
				print OUT "$total_reads\n";
			}else{
				my $mean_cov = int($total_reads/($bed_array[2]-$bed_array[1]));
				print OUT "$mean_cov\n";
			}
			$total_reads = 0;
			$bed_line = <BED>;
		}elsif((($bed_array[0] eq $pileup_array[0]) && ($pileup_array[1] <= $bed_array[1])) || ($chr_order{$pileup_array[0]} < $chr_order{$bed_array[0]})){
			# pileup position hasn't cought up to BED position yet, keep iterating
			$pileup_line = <STDIN>;
		}else{
			print "ERROR walking through files!! Problem spots:\n";
			print "BED: $bed_line\n";
			print "PILEUP: $pileup_line\n";
			last
		}
	}
}

close BED;
close OUT;






#!/usr/bin/perl

# This script will convert a sorted, binned bed file to wig
# Jared Evans evans.jared@mayo.edu

use strict;
use warnings;


if(scalar(@ARGV) != 5){
	die("USAGE: bed2wig.pl <input binned bed file> <exon bed file> <output wig file> <bin size> <coverage column>\n");
}

my $input = $ARGV[0];
my $exon_bed = $ARGV[1];
my $output = $ARGV[2];
my $bin = $ARGV[3];
my $coverage_col = $ARGV[4]-1;
	
open IN, "<$input" or die "opening $input\n";
open EXON_BED, "<$exon_bed" or die "opening $exon_bed\n";
open OUT, ">$output" or die "opening $output\n";

# prime pump
my $in_row = <IN>;
my $exon_row = <EXON_BED>;

while(defined $in_row || defined $exon_row){
	chomp $in_row if defined $in_row;
	chomp $exon_row if defined $exon_row;
	if(!defined $in_row){
		# no need to go further if reached end of BED
		last;
	}else{
		my @line = split("\t",$in_row);
		my $chr = $line[0];
		my $start = $line[1];
		my $stop = $line[2];
		my $coverage = $line[$coverage_col];
		if(defined $exon_row){
			# check to see if an exon header needs to be written
			my @exon_line = split("\t",$exon_row);
			if($chr eq $exon_line[0] and $start == $exon_line[1]){
				print OUT "fixedStep chrom=".$chr." start=".($start+1)." step=".$bin."\n";
				$exon_row = <EXON_BED>;
			}
		}
		print OUT $coverage."\n";
		$in_row = <IN>;
	}
}

close IN;
close EXON_BED;
close OUT;





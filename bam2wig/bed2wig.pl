#!/usr/bin/perl

# This script will convert a sorted, binned bed file to wig
# Jared Evans evans.jared@mayo.edu

use strict;
use warnings;


if(scalar(@ARGV) != 4){
	print "USAGE: bed2wig.pl <input binned bed file> <output wig file> <bin size> <coverage column>\n";
}else{
	my $input = $ARGV[0];
	my $output = $ARGV[1];
	my $bin = $ARGV[2];
	my $coverage_col = $ARGV[3]-1;
	
	open IN, "<$input" or die "opening $input\n";
	open OUT, ">$output" or die "opening $output\n";
	
	my $prev_start = 0;
	my $prev_chr = "";
	while(<IN>){
		my $row = $_;
		chomp $row;
		my @line = split("\t",$row);
		my $chr = $line[0];
		my $start = $line[1];
		my $stop = $line[2];
		my $coverage = $line[$coverage_col];
		$prev_chr = $chr if $prev_chr eq "";
		$prev_start = 0 if $prev_chr ne $chr; # handle chr transitions
		# print wig header for each exon
		if($start-$prev_start > $bin){
			print OUT "fixedStep chrom=".$chr." start=".($start+1)." step=".$bin."\n";
		}
		$prev_start = $start;
		$prev_chr = $chr;
		print OUT $coverage."\n";
		
		
	}
	
	close IN;
	close OUT;

	
}



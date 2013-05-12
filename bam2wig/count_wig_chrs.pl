#!/usr/bin/perl

# This script will count the total reads in each chr for a wig file
# This can be used as a QC method to ensure there isn't data missing from your wig
# Jared Evans evans.jared@mayo.edu

use strict;
use warnings;


if(scalar(@ARGV) != 1){
	print "USAGE: count_wig_chrs.pl <input wig file> \n";
}else{
	my $input = $ARGV[0];
	
	open IN, "<$input" or die "opening $input\n";
	
	my $wigname = $input;
	$wigname =~ s{.*/}{};      # removes path  
	#print $wigname;

	my $current_chr = "";
	my $current_chr_count = 0;
	while(<IN>){
		my $row = $_;
		chomp $row;
		my @line = split(" ",$row);
		my $reads = $line[0];
		if($reads eq "fixedStep") {
			my @chrom = split("=",$line[1]);
			if($chrom[1] ne $current_chr){
				print $current_chr."\t".$current_chr_count."\n" if $current_chr ne "";
				$current_chr = $chrom[1];
				$current_chr_count = 0;
			}
		}else{
			$current_chr_count+=$reads;
		}
	}

	print $current_chr."\t".$current_chr_count."\n" if $current_chr ne "";
	print "\n";
	close IN;

	
}



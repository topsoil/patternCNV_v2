#!/usr/bin/perl

# This script will split regions in a bed file into defined bin lengths
# When binning, if there are leftover lengths that are shorter than the bin size then they are ignored
# Jared Evans evans.jared@mayo.edu

use strict;
use warnings;


if(scalar(@ARGV) != 3){
	print "USAGE: bin_exons.pl <input bed file> <output BED file> <bin size>\n";
}else{
	my $input = $ARGV[0];
	my $output = $ARGV[1];
	my $bin = $ARGV[2];

	open IN, "<$input" or die "opening $input\n";
	open OUT, ">$output" or die "opening $output\n";
	
	
	while(<IN>){
		my $row = $_;
		chomp $row;
		my @line = split("\t",$row);
		my $chr = $line[0];
		my $start = $line[1];
		my $stop = $line[2];
		my $gene = "";
		$gene = "\t".$line[3] if defined $line[3];
		my $length = $stop-$start;
		# always print at least one bin for an exon
		#print OUT $chr."\t".$start."\t".($start+$bin).$gene."\n";
		#$start+=$bin;
		#$length-=$bin;
		while($length >= $bin){
			print OUT $chr."\t".$start."\t".($start+$bin).$gene."\n";
			$start+=$bin;
			$length-=$bin;
		}
	}
	
	close IN;
	close OUT;

}



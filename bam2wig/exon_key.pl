#!/usr/bin/perl

# Reformat and put finishing touches on exon key file
# Jared Evans evans.jared@mayo.edu

use strict;
use warnings;


if(scalar(@ARGV) != 3){
	print "USAGE: exon_key.pl <input exon bed, intersected with capture> <output exon key txt> <bin size>\n";
}else{
	my $input = $ARGV[0];
	my $output = $ARGV[1];
	my $bin_size = $ARGV[2];
	
	my $chr_col = 0;
	my $start_col = 1;
	my $stop_col = 2;
	my $gene_col = 3;
	my $intersect_col = 4;
	
	open IN, "<$input" or die "opening $input\n";
	open OUT, ">$output" or die "opening $output\n";
	
	# print header
	print OUT "Chr\tStart\tStop\tBin_Count\tGenes\tInCapture\n";
	while(<IN>){
		my $row = $_;
		chomp $row;
		my @line = split("\t",$row);
		my $chr = $line[$chr_col];
		my $start = $line[$start_col];
		my $stop = $line[$stop_col];
		my $genes = $line[$gene_col];
		my $incapture = $line[$intersect_col];
		
		#my $bin_count = int(($stop-$start)/$bin_size);
		
		#my $leftover = ($stop-($start+1))%$bin_size;
		my $leftover = ($stop-$start)%$bin_size;
		$stop=$stop-$leftover;
		my $bin_count = int(($stop-$start)/$bin_size);

		my @gene_list = split(";",$genes);
		my %uniq_genes = ();
		foreach my $gene (@gene_list){
			$uniq_genes{$gene} = 1;
		}
		my $uniq_gene_list = join(":",(sort keys %uniq_genes));
		$incapture = 1 if $incapture > 0;
		print OUT $chr."\t".($start+1)."\t".$stop."\t".$bin_count."\t".$uniq_gene_list."\t".$incapture."\n";
		
	}
	
	close IN;
	close OUT;

	
}



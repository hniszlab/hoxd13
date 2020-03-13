#!/usr/bin/perl 

use strict;
use warnings;

my %h;
my @l;
## column to use as index value
my $c=0;
my $num=0;
my @files=();

my $pv=shift;


foreach my $f(@ARGV){
	if($f =~ /(\d+)/){
		$num=$1;
	}else{
		$num++;
	}
	push(@files,$num);

	open IN,$f or die "Could not open file $f\n";
	while(<IN>){
		next if(/p_val/);
		chomp;
		@l=split();
		@{$h{$l[$c]}{$num}}=@l;
	}
	close IN;
}



my $val;
for my $gene(keys %h){
	print $gene;
	foreach my $f(sort {$a <=> $b} @files){
		$val=${$h{$gene}{$f}}[2];
		if($val){
			## adjusted p >= 0.05
			if(${$h{$gene}{$f}}[5] >= $pv){
				$val=0;
			}
		}
		$val||=0;

		## if adjusted p is greater than 0.05 then set to no change
		print "\t$val"
	}
	print "\n";
}




#!/usr/bin/perl

use strict;
use warnings;

sub round{
	my ($n,$d) = (@_);
	my $fac=10**$d;
	return (int($n*$fac)/$fac);
}


my $verbose=0;
$verbose=1 if($ARGV[4]);

open IN,$ARGV[0] or die "file1 not found\n";

my $ina;
my $inb;
my %in;
my $overlap=0;
my %out;

my $c1=0;
my $c2=0;

##define columns in file 1 and file 2 to combine (start counting from 1)
$c1=$ARGV[2]-1 if(defined $ARGV[2]);
$c2=$ARGV[3]-1 if(defined $ARGV[3]);

while(<IN>){
    chomp;
	next if(/^\s*$/);
	next if(/^#/);
	my @l=split();

	my @l2=split(",",$l[$c1]);
	foreach my $gene(@l2){
		$gene=$1 if($gene =~ /(\S+)_[ab]XXX/);
		next if($in{$gene}); ## no double counting

		$in{$gene} =1;
		$out{$gene}=$_;
		$ina++;
	}

}
close IN;

if($ARGV[5]){
	open OUT,">only_file2.txt" or die "Could not open file\n";
}
if($ARGV[5]){
	open OUT2,">file2_overlap_with_file1.txt" or die "Could not open file\n";
}


open IN,$ARGV[1] or die "file2 not found\n";
while(<IN>){
    chomp;
	next if(/^\s*$/);
	next if(/^#/);
	$inb++;
	my @l=split();
	if($in{$l[$c2]}){
		next if($in{$l[$c2]}>1); ## no double counting
		$overlap++; 
		print STDERR "o $_\n" if($verbose);
		$in{$l[$c2]}++;
		if($ARGV[5]){
			print OUT2 "$_\n";
		}
	}else{
		print STDERR "n $_\n" if($verbose);
		if($ARGV[5]){
			print OUT "$_\n";
		}

	}
}
close IN;
if($ARGV[5]){
close OUT;
}


print STDERR "#f1\t#f2\toverlap\tf1only\tf2only\tfrac_ov_1\tfrac_ov2\n";
print STDERR "$ina\t$inb\t$overlap\t",$ina-$overlap,"\t",$inb-$overlap,"\t",round($overlap/$ina,5),"\t",round($overlap/$inb,5),"\n";


if($ARGV[5]){
	open OUT,">overlap_files";
	for my $k(keys %out){
		print OUT "$out{$k}\n" if($in{$k} > 1);
	}
	close OUT;
}



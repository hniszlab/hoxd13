#!/usr/bin/env perl 

open IN,$ARGV[0] or die "File $ARGV[0] not found $!\n";

my $type='enhancer';

if($ARGV[1]){
	$type=$ARGV[1];
}

my $a=0;
while(<IN>){
	$a++;
	my @F=split();
	$F[1]++;
	print "$F[0]\t$F[3]_$a\t$type\t$F[1]\t$F[2]\t.\t.\t.\tID=$F[3]_$a\n";
}


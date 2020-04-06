#!/usr/bin/perl

open IN,$ARGV[0] or die "File not found\n";
open OUT,">$ARGV[1]" or die "Could not create $ARGV[1] file\n";
my $bc;
while(<IN>){
	@l=split();
	$bc='';
	if(/CB:Z:(\S+)/){
		$bc=$1;
	}
	$id="$l[0]:$bc";
	next if($seen{$id});
	$seen{$id}=1;

	print OUT "\@$id\n$l[9]\n+\n$l[10]\n";
}
close OUT;
close IN;

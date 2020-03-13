#!/usr/bin/perl 

open IN,$ARGV[0] or die "File not given\n";


my $cline=0;
my $offset=0;
while(<IN>){
	if(/start=(\d+)/){
		$offset=$1;
	}else{
		$cline++;
		while($cline < $offset){
			print "0\n";
			$cline++;
		}
		print;
	}
}
			


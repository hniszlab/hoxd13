#!/usr/bin/perl
use strict;
use warnings;


my %h;
## read in ids to check
foreach my $i(@ARGV){
	$h{$i}=1;
}

my $running=1;
while($running){
	$running=0;
	for my $k(keys %h){

		## grep for ids
		my $ret=`ps ax |grep $k`;
		my @l=split("\n",$ret);
		my $f=0;
		
		## check all grepped ids 
		foreach my $e(@l){
			my $id=$1 if($e =~ /^\s*(\d+)\s+/);
			$f=1 if($k == $id);
		}
		if($f == 0){
			delete $h{$k};
		}else{
			$running++;
		}
	}
	print STDERR "Still running $running jobs\r";
	sleep(10);
}
print STDERR "All jobs finished\n";


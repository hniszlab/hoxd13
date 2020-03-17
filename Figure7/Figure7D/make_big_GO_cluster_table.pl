#!/usr/bin/perl 

my @F=<GO_cluster*.csv>;
foreach my $f(@F){
	my $cl=$1 if($f =~ /(\d+)/);
	print "#cluster $cl\n";
	open IN,$f or die "Could not open file $f for reading\n";
	while(<IN>){
		print;
	}
	close IN;
}


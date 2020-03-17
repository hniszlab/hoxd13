#!/usr/bin/perl 

my @FILES=<*.csv>;

my %h;

my $cluster=0;

foreach my $f(@FILES){
	if($f=~ /(\d+)/){
		$cluster=$1;
	}
	open IN,$f or die "Could not open file $f for reading\n";
	while(<IN>){
		next if(/Enrichment/);
		@l=split(/\t/);
		$l[1] =~ s/ /_/g;
		next if($l[9] ne "x");
		$h{$l[1]}{$cluster}{'qv'}=$l[3];
		$h{$l[1]}{$cluster}{'pv'}=$l[2];
		$h{$l[1]}{$cluster}{'cs'}=$l[6];
		$h{$l[1]}{$cluster}{'cin'}=$l[8];
	}
	close IN;
}

print "#process";
for my $j(1..4){
for(my $i=1;$i<8;$i++){
	print "\t$i";
}
}
print "\n";
for my $k(keys %h){
	print $k;
	for(my $i=1;$i<8;$i++){
		if(exists $h{$k}{$i}){
			print "\t$h{$k}{$i}{'qv'}";
		}else{
			print "\t1";
		}
	}
	for(my $i=1;$i<8;$i++){
		if(exists $h{$k}{$i}){
			print "\t$h{$k}{$i}{'pv'}";
		}else{
			print "\t1";
		}
	}
	for(my $i=1;$i<8;$i++){
		if(exists $h{$k}{$i}){
			print "\t$h{$k}{$i}{'cs'}";
		}else{
			print "\t0";
		}
	}
	for(my $i=1;$i<8;$i++){
		if(exists $h{$k}{$i}){
			print "\t$h{$k}{$i}{'cin'}";
		}else{
			print "\t0";
		}
	}
	print "\n";
}


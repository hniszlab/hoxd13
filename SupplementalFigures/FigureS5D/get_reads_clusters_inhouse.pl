#!/usr/bin/perl 

use strict;
use warnings;

my $type="spdh";
$type=$ARGV[1] if($ARGV[1]);

my %hfc;
for(my $i=1;$i<=11;$i++){
	open my $fh,'>',"${type}_$i.sam" or die "Could not open filehandle for ${type}_$i.sam\n";
	$hfc{$i}=$fh;
}
## the ghost variable must be given or any other string 
## the filehandle as $h{c} must be in a block with { } ... 

#print { $h{'c'} } $_;

#close $h{'c'};
my %h;
my @l;
my $num=0;
my @FILES=<$ARGV[0]*.tsv>;

foreach my $f(@FILES){
	$num=0;
	if($f =~ /cluster_(\d+).tsv/){
		$num=$1;
	}

	open IN,$f or die "No file with barcodes given as $f\n";
	while(<IN>){
		chomp;
		@l=split();
		$h{$l[0]}=$num;
	}
	close IN;
}

my $p="/project/meissner_seq_data/10x_scRNAseq/2018-12-05_Hnisz/cr_count_mm10_mpimg_L15205-1_SPDH12-5-limb/outs";
if($type eq "spdh"){
}else{
	$p="/project/meissner_seq_data/10x_scRNAseq/2018-12-05_Hnisz/cr_count_mm10_mpimg_L15204-1_WT12-5-limb/outs";
}

my %bc;

my $lcount=0;
open IN,"samtools view $p/possorted_genome_bam.bam|" or die "Coult not open file possorted_genome_bam.bam\n";
my $file;
while(<IN>){
	$file='';
    if(/CB:Z:(\S+)/){
		if(not $h{$1}){
			$bc{$1}++;
		}else{
			$file=$h{$1};
			print { $hfc{$file} } $_;
		}
    }
}
close IN;

for(my $i=1;$i<=11;$i++){
	close $hfc{$i};
}
open OUT,">notbc_$type.csv";
for my $k(keys %bc){
	print OUT "$k\t$bc{$k}\n";
}
close OUT;

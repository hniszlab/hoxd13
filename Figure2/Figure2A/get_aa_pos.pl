#!/usr/bin/perl

# get_aa_pos.pl - Perl script
# Copyright (C) 2018 Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
use strict;
use warnings;

##        hydrophobic     uncharg pos   neg special
my @aa=qw(A V I L M F Y W S T N Q R H K D E C G P); ## U is Selenocystein

my ($id,$count,$protc);
my $seq;

open IN,$ARGV[0] or die "File not found\n";
while(<IN>){
    chomp;
    next if(/^\s*$/);
    if(/>(\S+)/){
        $id=$1;
    }else{
        $seq.=$_;
    }
}
close IN;
print "slen=",length($seq),"\n";

my @l=split("",$seq);
my %h;
for(my $i=0; $i <= $#l; $i++){
    push(@{$h{$l[$i]}},$i+1);
}            


## making R script now to load
$count++;
print "sl=list()\n";
print "sl[[$count]]=list()\n";

$protc++; ## if we have just one protein this is fine too
my $aac=0;
print "sl[[$count]][[$protc]]=list()\n";
foreach my $a(@aa){
    $aac++;
    ## if AA not occuring in sequence then set to 0
    if(not $h{$a}){
        print "sl[[$count]][[$protc]][[$aac]]=c(0)\n";
        ## else print the aa positions
    }else{
        print "sl[[$count]][[$protc]][[$aac]]=c(".join(",",@{$h{$a}}).")\n";
    }
}
print "aa<-c('".join("','",@aa)."')\n";

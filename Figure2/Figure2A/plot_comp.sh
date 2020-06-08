#!/bin/bash

if [[ -z $1 ]];then 
    echo 
    echo No fasta input file given
    echo 
    exit
fi

AA=A
if [[ $2 ]];then
    AA=$2
fi

perl get_aa_pos.pl $1 > tmp.R
Rscript amino_acid_comp_plot.R tmp.R $AA

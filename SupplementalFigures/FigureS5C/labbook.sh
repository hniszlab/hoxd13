#!/bin/bash
P1=../../Figure5/ScRNA-Seq_data_processing/marker_genes

j=0;
for i in {0..10};do 
        let j='i+1';echo cluster $j; 
		perl ../../src/a_in_b.pl $P1/mean_exprs_ctrl_clusterN_$j.tsv ../FigureS5B/cluster8_genes.tsv;done &>> overlaps_cluster9.txt

grep -P "^\d" overlaps_cluster9.txt  > overlaps_cluster9.nums
Rscript generate_figure_S5C.R


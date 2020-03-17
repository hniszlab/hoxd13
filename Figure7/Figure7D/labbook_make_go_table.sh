#!bin/bash

cut -f2,4 -d" " outmatrix_v4.csv |tail -n+2 > all_TFs


for i in {1..7};do 
	export CLUSTER=$i;perl -ane 'if($F[1] == $ENV{'CLUSTER'}){print "$F[0]\n";}' all_TFs > all_TFs_cluster_$i;done

cut -f1 -d" " all_TFs > all_TFs.genes

## now go to GORILLA and run GO analysis
## take all genes in data and the cluster specific ones. Save xls results
## open xls and save as csv.

## then run this one here
perl make_go_table.pl  > go_outtable.txt

## make big table for paper
perl make_big_GO_cluster_table.pl > big_GO_table.csv


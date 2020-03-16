#!/bin/bash 

for i in {1..22} X Y M;do
    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way/chr$i.phyloP100way.wigFix.gz .
done

## adjust and make nice pseudo random access files we can read in with fread from the data.table package
for i in {1..22};do
        perl transform_phyloP.pl chr${i}.phyloP100way.wigFix > chr${i}.phyloP100way.wigFix_ra
done


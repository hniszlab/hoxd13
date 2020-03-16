#!/bin/bash

## change the thread number here according to your needs
THREADS=1
## adjust he GENOME variable so that it points to the bowtie index file of mm9
## mm9 is here the prefix of the index file. E.g. mm9.1.ebwt would be the first index file
GENOME=mm9

## SRR3498935 is the hoxd13 ChIP
## SRR3498935 is the hoxd13 INPUT 

## mapping
for num in 35 39;do
	FILE=SRR34989${num}
	BP='-n 2 -e 70 -m 1 -k 1 --best -l 200'
	bowtie --chunkmbs 200 --threads=$THREADS --sam $BP $GENOME <(gunzip -dc ${FILE}.fastq.gz) $FILE.sam
done

## running macs14
macs14 -t SRR3498935.sam -c SRR3498939.sam -g mm9 --pvalue 10 -f SAM -w -S --space=50 --keep-dup=auto

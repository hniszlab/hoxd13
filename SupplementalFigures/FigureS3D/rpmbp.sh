#!/bin/bash

if [[ -z $2 ]];then
	echo
	echo USAGE: bash rpmbp.sh gff  bam  outfilename  bin_number[def=1000] extension[def=0]
	echo
	exit
fi


## path to BradnerPipeline. Make sure it points to the right location
P=~/a/tools/pipeline


CWD=$(pwd -LP)
ANNO=$(readlink -f $1)
bn=$(basename $ANNO)

if [[ $ANNO =~ \.bed$ ]];then
	bn=$(basename $ANNO .bed)
	if [[ ! -f $bn.gff ]];then
		perl bed2gff.pl $ANNO > $bn.gff
	fi
	ANNO=$(readlink -f $bn.gff)
fi


BAM=$2
if [[ ! -f $BAM ]];then
	echo bam file not found
	exit
fi

DATE=$(date +%Y%m%d%k%M%S)
OUT=$bn.rpmbp_$DATE
if [[ $3 ]];then
	OUT=$3
	if [[ -f $OUT ]];then 
		OUT="${3}_$DATE"
		echo outfile is $OUT since file $3 exists already
	fi
fi


BINS=1000
if [[ $4 ]];then
	BINS=$4
fi
EXT=0

if [[ $5 ]];then
	EXT=$5
fi

echo bamToGFF_turbo.py -b $BAM -i $ANNO -o $CWD/$OUT -s both -e $EXT -m $BINS -r
if [[ ! -f "$CWD/$OUT" ]];then
	python2 $P/bamToGFF_turbo.py -b $BAM -i $ANNO -o $CWD/$OUT -s both -e $EXT -m $BINS -r
else
	echo File $CWD/$OUT exists
	exit
fi

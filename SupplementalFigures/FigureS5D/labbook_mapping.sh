#!/bin/bash

STARINDEX=/star/mm9

## barcodes per cluster are generated in the Figure5/ScRNA-Seq folder
BC=../../Figure5/ScRNA-Seq_data_processing/barcodes/

## getting reads first 
for tt in spdh ctrl;do 
	echo $tt
	perl get_reads_clusters.pl $BC/barcodes_$tt $tt
done

## make fastq files ready for mapping
for tt in spdh ctrl;do 
	for i in {1..11};do 
		echo $i
		perl sam_to_fastq.pl ${tt}_${i}.sam ${tt}_${i}.fastq
	done
done


## map to mm9 
for i in {1..11};do 
	THREADS=50
	echo $i
	STAR --genomeDir $STARINDEX  --readFilesIn ctrl_$i.fastq --runThreadN $THREADS --outFileNamePrefix ctrl_${i}_mm9.sam
done

for i in {1..11};do 
	THREADS=50
	echo $i 
	STAR --genomeDir $STARINDEX  --readFilesIn spdh_$i.fastq --runThreadN $THREADS --outFileNamePrefix spdh_${i}_mm9.sam
done


## merge bam files
bamtools merge -list spdh_bam.list -out spdh_merged.bam
bamtools merge -list ctrl_bam.list -out ctrl_merged.bam

## need to build an index first
for tt in ctrl spdh;do
	THREADS=50
	FILE=${tt}_merged.bam
	bn=$(basename $FILE .bam)
	samtools view --threads $THREADS -bS $FILE | samtools sort --threads $THREADS - > ${bn}_sorted.bam
	samtools index -@ $THREADS ${bn}_sorted.bam ${bn}_sorted.bai
done


BIN=$HOME/bin/deeptools
THREADS=8

## normalize, for paper we use BPM files
for NORM in CPM BPM; do
	for tt in ctrl spdh;do 
		for i in merged;do
			ID=${tt}_${i}_${NORM}
			if [[ ! -f cluster_${ID}.bw ]];then
				$BIN/bamCoverage -b ${tt}_${i}_sorted.bam -o cluster_${ID}.bw --normalizeUsing $NORM -p $THREADS
			fi
		done
	done
done

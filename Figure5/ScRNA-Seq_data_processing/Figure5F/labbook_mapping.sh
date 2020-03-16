#!/bin/bash

## getting reads first 
for tt in spdh ctrl;do 
	echo $tt
	perl get_reads_clusters.pl barcodes_$tt $tt
done

## make fastq files ready for mapping
for tt in spdh ctrl;do 
for i in {1..11};do 
	echo $i
	perl sam_to_fastq.pl ${tt}_${i}.sam ${tt}_${i}.fastq
done
done


for i in {1..11};do echo $i;STAR --genomeDir /project/bioinf_meissner/references/mm9  --readFilesIn ctrl_$i.fastq --runThreadN 50 --outFileNamePrefix ctrl_${i}_mm9.sam;done
for i in {1..11};do echo $i;STAR --genomeDir /project/bioinf_meissner/references/mm9  --readFilesIn ~/a/scRNAseq_limb/get_reads_per_cluster/spdh_$i.fastq --runThreadN 50 --outFileNamePrefix spdh_${i}_mm9.sam;done


## merge bam files
time bamtools merge -list spdh_bam.list -out spdh_merged.bam
time bamtools merge -list ctrl_bam.list -out ctrl_merged.bam

## need to build an index first

for tt in ctrl spdh;do
THREADS=50
FILE=${tt}_merged.bam
bn=$(basename $FILE .bam)

samtools view --threads $THREADS -bS $FILE | samtools sort --threads $THREADS - > ${bn}_sorted.bam
samtools index -@ $THREADS ${bn}_sorted.bam ${bn}_sorted.bai

done


## make tracks 
LOG=/project/hnisz_lab_analyses/logs
PP=/project/hnisz_lab_analyses/scRNAseq_limb/get_reads_per_cluster
BIN=/project/hnisz_lab_analyses/tools/bin/deeptools


for NORM in CPM BPM; do
for tt in ctrl spdh;do for i in merged;do
    ID=${tt}_${i}_${NORM}
    if [[ ! -f cluster_${ID}.bw ]];then
        mxqsub --stdout=$LOG/stdout_$ID.log --stderr=$LOG/stderr_$ID.log --threads=8 --memory=20000 --runtime=480  $BIN/bamCoverage -b $PP/${tt}_${i}_sorted.bam -o $PP/cluster_${ID}.bw --normalizeUsing $NORM -p 8
    fi
    done
done
done

cp *merged*.bw /project/smdn1801/httpd-2.4.34/htdocs/track_hubs/mm9_rpmbp/

smdnp
/project/smdn1801/httpd-2.4.34/htdocs/track_hubs/mm9_rpmbp

for i in `ls *merged_BPM.bw`;do echo $i|perl -ane 'if(/(\d+)/){$num=$1;};chomp;@l=split("_BPM");print "$_\t$l[0]\tWT_vs_SPDH_cluster_merged\n";'; done >> nfile_merged

## make merged version
perl make_trackdb_overlay_full_vis_and_autoscale.pl nfile_merged 1 > trackDb_cluster_merged.txt

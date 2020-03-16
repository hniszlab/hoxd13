#!/bin/bash


DPATH=../../Fdata/FigureS3D/

bamWT=$DPATH/3749_3750_v3_gg3_merged_sorted.bam
bam7A=$DPATH/3746_GCCAAT_R1_001_v3_gg3_merged_sorted.bam

bedwt=../FigureS3D/rpmbp_comp_wt12_vs_7A/wt12_vs_7A.intersect_tag_fused_sorted.bed

## this is for the union of all peaks
bash meta_plots.sh $bedwt $bamWT wt_vs_7A_wt12 &
bash meta_plots.sh $bedwt $bam7A wt_vs_7A_7A &

## bam files
bamWTpub=$DPATH/HOXD13r12_v3_gg3_merged_sorted.bam 
bamQR=$DPATH/HOXD13QRr12_v3_gg3_merged_sorted.bam

## bed files
bed_WT_QR=../FigureS3D/rpmbp_comp_D13r12_vs_QR/D13r12_vs_QR.intersect_tag_fused_sorted.bed


## run meta analyis now
bash meta_plots.sh $bed_WT_QR $bamWTpub  R1R2_vs_QR_wtR1R2 &
bash meta_plots.sh $bed_WT_QR $bamQR        R1R2_vs_QR_QR &


for COLOR in BuRd;do
	Rscript merge_plots.R rpmbp_wt_vs_7A_wt12 rpmbp_wt_vs_7A_7A wt12_vs_7A 3000 $COLOR
	Rscript merge_plots.R rpmbp_R1R2_vs_QR_wtR1R2 rpmbp_R1R2_vs_QR_QR D13r12_vs_QR 3000 $COLOR
done



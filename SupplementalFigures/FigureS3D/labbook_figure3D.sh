#!/bin/bash

DPATH=../../Fdata/FigureS3D/

#PV=5
python3 compare_bed_files_rpmbp.py -a $DPATH/3749_3750_v3_gg3_merged_sorted_5_peaks.bed \
	-A $DPATH/3746_GCCAAT_R1_001_v3_gg3_merged_sorted_5_peaks.bed \
	-t wt12 -T 7A \
	-b $DPATH/3749_3750_v3_gg3_merged_sorted.bam \
	-B $DPATH/3746_GCCAAT_R1_001_v3_gg3_merged_sorted.bam


#PV=10
python3 compare_bed_files_rpmbp.py -a $DPATH/HOXD13r12_v3_gg3_merged_sorted_5_peaks.bed \
	-A $DPATH/HOXD13QRr12_v3_gg3_merged_sorted_5_peaks.bed \
	-t D13r12 -T QR \
	-b $DPATH/HOXD13r12_v3_gg3_merged_sorted.bam \
	-B $DPATH/HOXD13QRr12_v3_gg3_merged_sorted.bam


## making paper figures, Note. The 
cd rpmbp_comp_wt12_vs_7A
Rscript ../make_raster_plot.R rpmbp_wt12_vs_7A.intersect_tag_fused_sorted.bed_wt12 rpmbp_wt12_vs_7A.intersect_tag_fused_sorted.bed_7A HOXD13_merged SPDH
cd -
cd rpmbp_comp_D13r12_vs_QR
Rscript make_raster_plot.R rpmbp_D13r12_vs_QR.intersect_tag_fused_sorted.bed_D13r12 rpmbp_D13r12_vs_QR.intersect_tag_fused_sorted.bed_QR HOXD13_merged HOXD13_QR
cd -



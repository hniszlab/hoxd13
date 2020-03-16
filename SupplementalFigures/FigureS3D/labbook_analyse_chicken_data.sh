#!/bin/bash
## mapping

## these are the files from daniel

GENOME=gg3

for j in *.gz;do 
	i=$(basename $j fastq.gz); 
	mxqsub --stdout=$LOG/stdout_$i.log --stderr=$LOG/stderr_$i.log --threads=8 --memory=6000 --runtime=600  $BIN/runbowtie.sh $GENOME $j 1 0 $GENOME 1
done


## these are the files from the publication
cd 
for i in `cat ~/filesibra`;do 
	mxqsub --stdout=$LOG/stdout_$i.log --stderr=$LOG/stderr_$i.log --threads=8 --memory=6000 --runtime=600  $BIN/runbowtie.sh $GENOME $i 1 0 $GENOME 1
done

## waiting until this is finished


for i in `ls 3746_*${GENOME}*.sam`;do sam_to_bigwig.sh $i 50 1;done
for i in `ls 3749_*${GENOME}*.sam`;do sam_to_bigwig.sh $i 50 1;done
for i in `ls 3750_*${GENOME}*.sam`;do sam_to_bigwig.sh $i 50 1;done

## we need to merge lanes now to a single file 



GS=1046932099
PV=5
for j in 2 3;do

bamtools merge -in 3750_AGTTCC_L001_R1_001_v${j}_${GENOME}_sorted.bam -in 3750_AGTTCC_L002_R1_001_v${j}_${GENOME}_sorted.bam -out 3750_AGTTCC_R1_001_v${j}_${GENOME}_merged.bam
sam_to_bigwig.sh 3750_AGTTCC_R1_001_v${j}_${GENOME}_merged.bam 50 1
bamtools merge -in 3749_AGTCAA_L001_R1_001_v${j}_${GENOME}_sorted.bam -in 3749_AGTCAA_L002_R1_001_v${j}_${GENOME}_sorted.bam -out 3749_AGTCAA_R1_001_v${j}_${GENOME}_merged.bam
sam_to_bigwig.sh 3749_AGTCAA_R1_001_v${j}_${GENOME}_merged.bam 50 1
bamtools merge -in 3746_GCCAAT_L001_R1_001_v${j}_${GENOME}_sorted.bam -in 3746_GCCAAT_L002_R1_001_v${j}_${GENOME}_sorted.bam -out 3746_GCCAAT_R1_001_v${j}_${GENOME}_merged.bam
sam_to_bigwig.sh 3746_GCCAAT_R1_001_v${j}_${GENOME}_merged.bam 50 1
done
for j in 2 3;do
run_macs_single.sh 3750_AGTTCC_R1_001_v${j}_${GENOME}_merged_sorted NONE $GS $PV NONE BAM &
run_macs_single.sh 3749_AGTCAA_R1_001_v${j}_${GENOME}_merged_sorted NONE $GS $PV NONE BAM &
run_macs_single.sh 3746_GCCAAT_R1_001_v${j}_${GENOME}_merged_sorted NONE $GS $PV NONE BAM &
done

## make sorted bam files
bash make_ibrafiles.sh $GENOME

## run macs
cd ~/a/mapped/
bash ibramacs.sh $GENOME

## edit bed files
cd ~/a/macs14
## the second replicate seems to be lost totally
python3 add_id.py 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks.bed 7A Rep2 > 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed
python3 add_id.py 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks.bed WT Rep1 > 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed
python3 add_id.py 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks.bed WT Rep2 > 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed
python3 add_id.py SRR771381_${GENOME}_sorted_10_peaks.bed QK Rep1 > SRR771381_${GENOME}_sorted_10_peaks_QK_Rep1.bed
python3 add_id.py SRR771382_${GENOME}_sorted_10_peaks.bed QK Rep2 > SRR771382_${GENOME}_sorted_10_peaks_QK_Rep2.bed
python3 add_id.py SRR771383_${GENOME}_sorted_10_peaks.bed QR Rep1 > SRR771383_${GENOME}_sorted_10_peaks_QR_Rep1.bed

A=$HOME/a/macs14
M=$HOME/a/mapped
GS=1046932099
PV=10
GENOME=gg3

## wt
bed0_1=3749_AGTCAA_R1_001_v3_${GENOME}_merged_sorted_5_peaks.bed
bed0_2=3750_AGTTCC_R1_001_v3_${GENOME}_merged_sorted_5_peaks.bed
bam0_1=3749_AGTCAA_R1_001_v3_${GENOME}_merged_sorted.bam
bam0_2=3750_AGTTCC_R1_001_v3_${GENOME}_merged_sorted.bam

## 7A mutant
bed1=3746_GCCAAT_R1_001_v3_${GENOME}_merged_sorted_5_peaks.bed
bam1=3746_GCCAAT_R1_001_v3_${GENOME}_merged_sorted.bam

## now we take the official hoxd13 wildtypes 
bed2_1=SRR771378_${GENOME}_sorted_${PV}_peaks.bed
bed2_2=SRR771379_${GENOME}_sorted_${PV}_peaks.bed
bam2_1=SRR771378_${GENOME}_sorted.bam
bam2_2=SRR771379_${GENOME}_sorted.bam

#HOXD13_Q->K mutant
bed3_1=SRR771381_${GENOME}_sorted_${PV}_peaks_QK_Rep1.bed
bed3_2=SRR771382_${GENOME}_sorted_${PV}_peaks_QK_Rep2.bed
bam3_1=SRR771381_${GENOME}_sorted.bam
bam3_2=SRR771382_${GENOME}_sorted.bam

#HOXD13_Q->R mutant
bed4_1=SRR771383_${GENOME}_sorted_{$PV}_peaks.bed
bed4_2=SRR771384_${GENOME}_sorted_{$PV}_peaks.bed

#bed3_2=SRR771382_${GENOME}_sorted_10_peaks_QK_Rep2.bed
bam4_1=SRR771383_${GENOME}_sorted.bam
bam4_2=SRR771384_${GENOME}_sorted.bam
#bam3_2=SRR771382_${GENOME}_sorted.bam


# merge wt samples 
bamtools merge -in $bam0_1 -in $bam0_2 -out 3749_3750_v3_gg3_merged.bam
sam_to_bigwig.sh 3749_3750_v3_gg3_merged.bam 50 1
run_macs_single.sh 3749_3750_v3_gg3_merged_sorted NONE $GS $PV NONE BAM &

#compare_bed_files_rpmbp.py $A/3749_3750_v3_gg3_merged_sorted.bed $A/$bed1 wt1 7A $M/3749_3750_v3_gg3_merged_sorted.bam $M/$bam1

PV=10
GENOME=gg3
## wt from publication
bamtools merge -in $bam2_1 -in $bam2_2 -out HOXD13r12_v3_gg3_merged.bam
sam_to_bigwig.sh HOXD13r12_v3_gg3_merged.bam 50 1
run_macs_single.sh HOXD13r12_v3_gg3_merged_sorted NONE $GS $PV NONE BAM 

## QK from publication
bamtools merge -in $bam3_1 -in $bam3_2 -out HOXD13QKr12_v3_gg3_merged.bam
sam_to_bigwig.sh HOXD13QKr12_v3_gg3_merged.bam 50 1
run_macs_single.sh HOXD13QKr12_v3_gg3_merged_sorted NONE $GS $PV NONE BAM &

## QR from publication
bamtools merge -in $bam4_1 -in $bam4_2 -out HOXD13QRr12_v3_gg3_merged.bam
sam_to_bigwig.sh HOXD13QRr12_v3_gg3_merged.bam 50 1
run_macs_single.sh HOXD13QRr12_v3_gg3_merged_sorted NONE $GS $PV NONE BAM &

## maybe we should use the input controls as well at least for the publication ones
## the the second input for the QR is missing so not using for all three of them is then more fair also likely more dirty
## we could use the IDR tool actually to see what the values look like then as beeing reproducibly called with a pvalue
PV=10
SUFP=v3_gg3_merged_sorted_${PV}_peaks.bed
SUFB=v3_gg3_merged_sorted.bam
## now compare them 
compare_bed_files_rpmbp.py -a $A/3749_3750_$SUFP -A $A/$bed1 -t wt12 -T 7A -b $M/3749_3750_$SUFB -B $M/$bam1
compare_bed_files_rpmbp.py -a $A/$bed0_1 -A $A/$bed0_2 -t wt1 -T wt2 -b $M/$bam0_1 -B $M/$bam0_2

## new run ## should add pvalue to it so we dont mix it up at least
compare_bed_files_rpmbp.py -a $A/$bed0_1 -A $A/$bed0_2 -t wt1 -T wt2 -b $M/$bam0_1 -B $M/$bam0_2
## looks good


## wt published vs QR
## rerun those
## run macs on merged files too!!!! 


compare_bed_files_rpmbp.py -a $A/HOXD13r12_$SUFP -A $A/HOXD13QKr12_$SUFP -t D13r12 -T QK -b $M/HOXD13r12_$SUFB -B $M/HOXD13QKr12_$SUFB
compare_bed_files_rpmbp.py -a $A/HOXD13r12_$SUFP -A $A/HOXD13QRr12_$SUFP -t D13r12 -T QR -b $M/HOXD13r12_$SUFB -B $M/HOXD13QRr12_$SUFB

## redo, will give new outfile name as HOXD13r1_vs_HOXD13r2
compare_bed_files_rpmbp.py -a $A/$bed2_1 -A $A/$bed2_2 -t HOXD13r1 -T HOXD13r2 -b $M/$bam2_1 -B $M/$bam2_2

##

### make plots for paper
cd rpmbp_comp_wt1_wt2
Rscript ../plot.R rpmbp_wt1_vs_wt2.intersect_tag_fused_sorted.bed_wt1 rpmbp_wt1_vs_wt2.intersect_tag_fused_sorted.bed_wt2 HOXD13_Rep1 HOXD13_Rep2
cd ..
cd rpmbp_comp_wt12_7A
Rscript ../plot.R rpmbp_wt12_vs_7A.intersect_tag_fused_sorted.bed_wt12 rpmbp_wt12_vs_7A.intersect_tag_fused_sorted.bed_7A HOXD13_merged SPDH
cd ..

cd rpmbp_comp_HOXD13r1_HOXD13r2
Rscript ../plot.R rpmbp_HOXD13r1_vs_HOXD13r2.intersect_tag_fused_sorted.bed_HOXD13r1 rpmbp_HOXD13r1_vs_HOXD13r2.intersect_tag_fused_sorted.bed_HOXD13r2 HOXD13p_Rep1 HOXD13p_Rep2
cd ..
cd rpmbp_comp_D13r12_QR
Rscript ../plot.R rpmbp_D13r12_vs_QR.intersect_tag_fused_sorted.bed_D13r12 rpmbp_D13r12_vs_QR.intersect_tag_fused_sorted.bed_QR HOXD13_merged HOXD13_QR
cd ..
cd rpmbp_comp_D13r12_QK
Rscript ../plot.R rpmbp_D13r12_vs_QK.intersect_tag_fused_sorted.bed_D13r12 rpmbp_D13r12_vs_QK.intersect_tag_fused_sorted.bed_QK HOXD13_merged HOXD13_QK
cd ..


if [[ $abc ]];then
	compare_bed_files_rpmbp.py $A/$bed0_1 $A/$bed0_2 wt1 wt2 $M/$bam0_1 $M/$bam0_2
	compare_bed_files_rpmbp.py $A/$bed0_1 $A/$bed1 wt1 7A $M/$bam0_1 $M/$bam1
	compare_bed_files_rpmbp.py $A/$bed0_2 $A/$bed1 wt2 7A $M/$bam0_2 $M/$bam1
fi

#compare_bed_files_rpmbp.py bedfile1 bedfile2 tag1 tag2 bamfile1 bamfile2
compare_bed_files_rpmbp.py $A/$bed1 $A/$bed2_1 7A HOXD13r1 $M/$bam1 $M/$bam2_1
compare_bed_files_rpmbp.py $A/$bed1 $A/$bed2_2 7A HOXD13r2 $M/$bam1 $M/$bam2_2

## compare old wt with official
compare_bed_files_rpmbp.py $A/$bed0_1 $A/$bed2_1 wt1 HOXD13r1 $M/$bam0_1 $M/$bam2_1 &
compare_bed_files_rpmbp.py $A/$bed0_1 $A/$bed2_2 wt1 HOXD13r2 $M/$bam0_1 $M/$bam2_2 &
compare_bed_files_rpmbp.py $A/$bed0_2 $A/$bed2_1 wt2 HOXD13r1 $M/$bam0_2 $M/$bam2_1 &
compare_bed_files_rpmbp.py $A/$bed0_2 $A/$bed2_2 wt2 HOXD13r2 $M/$bam0_2 $M/$bam2_2 &

## the public chicken replicates 
compare_bed_files_rpmbp.py $A/$bed2_1 $A/$bed2_2 HOXD13r1 HOXD13r2 $M/$bam2_1 $M/$bam2_2


compare_bed_files_rpmbp.py $A/$bed2_1 $A/$bed3_1 HOXD13r1 QtoK1 $M/$bam2_1 $M/$bam3_1 &
compare_bed_files_rpmbp.py $A/$bed2_1 $A/$bed3_2 HOXD13r1 QtoK2 $M/$bam2_1 $M/$bam3_2 &
compare_bed_files_rpmbp.py $A/$bed2_2 $A/$bed3_1 HOXD13r2 QtoK1 $M/$bam2_2 $M/$bam3_1 &
compare_bed_files_rpmbp.py $A/$bed2_2 $A/$bed3_2 HOXD13r2 QtoK2 $M/$bam2_2 $M/$bam3_2 &


# compare_bed_files_rpmbp.py $A/$bed1 $A/$bed2_1 7A HOXD13r1 $M/$bam1 $M/$bam2_1


#python3 add_id.py SRR771384_${GENOME}_sorted_10_peaks.bed QR_Rep2 > SRR771384_${GENOME}_sorted_10_peaks_QR_Rep2.bed
python3 add_id.py SRR771378_${GENOME}_sorted_10_peaks.bed HOXD13 Rep1 > SRR771378_${GENOME}_sorted_10_peaks_HOXD13_Rep1.bed
python3 add_id.py SRR771379_${GENOME}_sorted_10_peaks.bed HOXD13 Rep2 > SRR771379_${GENOME}_sorted_10_peaks_HOXD13_Rep2.bed

## wt 1 vs wt2
cat 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed > wt_R1_vs_wt_R2.bed


## mutants vs each other
cat 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed SRR771381_${GENOME}_sorted_10_peaks_QK_Rep1.bed > ALA7_QK_R1.bed
cat 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed SRR771382_${GENOME}_sorted_10_peaks_QK_Rep2.bed > ALA7_QK_R2.bed
cat 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed SRR771383_${GENOME}_sorted_10_peaks_QR_Rep1.bed > ALA7_QR_Rep1.bed
#cat 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed SRR771384_${GENOME}_sorted_10_peaks_QR_Rep2.bed > ALA7_QR_Rep2.bed

## all wt vs 7A
cat 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed > wt_R1_vs_ALA7.bed
cat 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed > wt_R2_vs_ALA7.bed
cat SRR771378_${GENOME}_sorted_10_peaks_HOXD13_Rep1.bed 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed > hoxd13_R1_vs_ALA7.bed
cat SRR771379_${GENOME}_sorted_10_peaks_HOXD13_Rep2.bed 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed > hoxd13_R2_vs_ALA7.bed

## all wt vs QK
VS=SRR771381_${GENOME}_sorted_10_peaks_QK_Rep1.bed
cat 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed $VS > wt_R1_vs_QK_R1.bed
cat 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed $VS > wt_R2_vs_QK_R1.bed

VS=SRR771382_${GENOME}_sorted_10_peaks_QK_Rep2.bed
cat 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed $VS > wt_R1_vs_QK_R2.bed
cat 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed $VS > wt_R2_vs_QK_R2.bed

## wt vs QR
VS=SRR771383_${GENOME}_sorted_10_peaks_QR_Rep1.bed
cat 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed $VS > wt_R1_vs_QR_R1.bed
cat 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed $VS > wt_R2_vs_QR_R1.bed
cat SRR771378_${GENOME}_sorted_10_peaks_HOXD13_Rep1.bed 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed > hoxd13_R1_vs_ALA7.bed
cat SRR771379_${GENOME}_sorted_10_peaks_HOXD13_Rep2.bed 3746_GCCAAT_R1_001_v3_${GENOME}_sorted_merged_5_peaks_7A.bed > hoxd13_R2_vs_ALA7.bed


## QR vs ALA7
cat SRR771383_${GENOME}_sorted_10_peaks_QR_Rep1.bed 


## all WT versus each other
cat SRR771378_${GENOME}_sorted_10_peaks_HOXD13_Rep1.bed 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed > hoxd13_R1_vs_wt_R1.bed
cat SRR771379_${GENOME}_sorted_10_peaks_HOXD13_Rep2.bed 3749_AGTCAA_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep1.bed > hoxd13_R2_vs_wt_R1.bed
cat SRR771378_${GENOME}_sorted_10_peaks_HOXD13_Rep1.bed 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed > hoxd13_R1_vs_wt_R2.bed
cat SRR771379_${GENOME}_sorted_10_peaks_HOXD13_Rep2.bed 3750_AGTTCC_R1_001_v3_${GENOME}_sorted_merged_5_peaks_WT_Rep2.bed > hoxd13_R2_vs_wt_R2.bed

## all WT vs 

## make gff file
bed2gff.py wt_R1_vs_ALA7.bed > wt_R1_vs_ALA7.gff
bed2gff.py wt_R2_vs_ALA7.bed > wt_R2_vs_ALA7.gff
bed2gff.py wt_R1_vs_QK_R1.bed > wt_R1_vs_QK_R1.gff 
bed2gff.py wt_R2_vs_QK_R1.bed > wt_R2_vs_QK_R1.gff
bed2gff.py wt_R1_vs_QK_R1.bed > wt_R1_vs_QK_R2.gff 
bed2gff.py wt_R2_vs_QK_R2.bed > wt_R2_vs_QK_R2.gff
bed2gff.py hoxd13_R1_vs_wt_R1.bed > hoxd13_R1_vs_wt_R1.gff
bed2gff.py hoxd13_R1_vs_wt_R2.bed > hoxd13_R1_vs_wt_R2.gff
bed2gff.py hoxd13_R2_vs_wt_R1.bed > hoxd13_R2_vs_wt_R1.gff
bed2gff.py hoxd13_R2_vs_wt_R2.bed > hoxd13_R2_vs_wt_R2.gff
bed2gff.py wt_R1_vs_wt_R2.bed > wt_R1_vs_wt_R2.gff
bed2gff.py ALA7_QK_R1.bed >ALA7_QK_R1.gff 
bed2gff.py ALA7_QK_R2.bed >ALA7_QK_R2.gff 

## run rose
DIR=$(pwd -LP)
TYPE=se
EXT=200
DIR=$(pwd -LP)
BINS=50
A=/project/hnisz_lab_analyses/rose/h3k27acWT_vs_h3k27ac_input_pv_1e-8_tss_2_5kb_se
AN=/project/hnisz_lab_analyses
MACS=$AN/macs14

## HOX13 acetyl with windows
cd $AN/tools/pipeline

#for B in *bam;do
#python bamToGFF_turbo.py -b $AN/mapped/$B -i $DIR/$bn.gff -o $DIR/rpmbp_${bn}_centerPeak_${Bbn}_${BINS}bins -s both -e $EXT -m $BINS -r &

## run it for all necessary combinations
cat ~/s/daniel/files_for_rpmbp |grep -v \\#|perl -ane 'if(/\S/){print "cd \$AN/tools/pipeline\npython bamToGFF_turbo.py -b \$AN/mapped/$F[1] -i \$DIR/$F[0] -o \$DIR/rpmbp_${F[2]}_centerPeak_\${Bbn}_\${BINS}bins -s both -e \$EXT -m \$BINS -r \&\ncd \$DIR\n";}'

cd $DIR
for i in *centerPeak__50bins;do cut -f3-27 $i|tail -n+2 > $i.csv; done
## now the plots
## make plots

#!/bin/bash 
AN=/project/hnisz_lab_analyses/
W=(pwd -LP)

## all promoters bed file 
ALLPROM=/project/hnisz_lab_storage/annotation/mm9_refseq_20180924_tss_refseq_symbol_w1000nt_active_DI_t-1.bed

## get tads in the limb
TAD=$AN/annotation/mouse_limb_E11.5_mm9_all_mapq30_KR_250kb_all_domains_TADs_corrected.bed

## more TADs
perl -ane 'if($F[2]-$F[1] < 1e6){print;}' $TAD > tads_1mb.bed
perl -ane 'if($F[2]-$F[1] < 5e5){print;}' $TAD > tads_0.5mb.bed
perl -ane 'if($F[2]-$F[1] < 1e5){print;}' $TAD > tads_0.1mb.bed

wc -l tads_*mb.bed
#    67 tads_0.1mb.bed
#  2862 tads_0.5mb.bed
#  4028 tads_1mb.bed

TAD=tads_1mb.bed

### we decided to use the single end mapping only and then of course only read 1 ?
#HOXPE=$HOME/a/macs14/SRR3498935_peaks.bed

## HOXD13 peaks are E11.5 
HOX=$HOME/a/macs14/SRR3498935_se_10_peaks.bed

## it were all cluster4 dysregulated genes, the naming here is misleading!!! 
#C4=cluster_4_up_in_spdh.bed
C4=cluster_4_dysregulated.bed
## get_gene_models back from bed file
#perl get_active_genes_gene_models_back.pl cluster_4_up_in_spdh.txt > $C4
perl get_active_genes_gene_models_back.pl  ../diffex_genes_cluster_4.tsv > $C4


## now intersect with TADs
intersectBed -a $TAD -b $C4 -wo -F 1.0 > cluster4inTAD
intersectBed -a $TAD -b $HOX -wo -F 1.0 > HoxinTAD

## control genes maybe
CTRLG=/project/hnisz_lab_storage/annotation/mm9_refseq_20180924_active_DI_t-1.bed

## get genes that are not in cluster4
a_in_b.pl  $C4 $CTRLG 4 4 0 1 

UNG=NM_not_cluster4.bed
grep -v NR_ only_file2.txt > NM_not_cluster4.bed

## take out gene symbols of our VP genes in cluster4, more correct since now we use TADs not based on TADs in the orig set
perl exclude.pl cluster_4_up_in_spdh.bed NM_not_cluster4.bed > NM_not_cluster4_all.bed
## what happens if we take all genes coding and noncoding
perl exclude.pl cluster_4_up_in_spdh.bed only_file2.txt > not_cluster4_all.bed

intersectBed -a $TAD -b $UNG -wo -F 1.0 > CtrlginTAD

#
### how many tads with genes
#perl fuse.pl $UNG CtrlginTAD HoxinTAD g |wc -l
##3567
### how many tads with genes and hoxd peak
#perl fuse.pl $UNG CtrlginTAD HoxinTAD g|perl -ane 'print if($F[1] >0);'|wc -l
##2760
#
### how many hox in tad
#perl fuse.pl $UNG CtrlginTAD HoxinTAD h |wc -l
#3369
#
### how many genes in thos TADS
#perl fuse.pl $UNG CtrlginTAD HoxinTAD g|perl -ane 'if($F[2] >0){$sum+=$F[2];};END{print "$sum\n";};'
#24650
#
### how many hoxd13 peaks in tads with $UNG genes
#perl fuse.pl $UNG CtrlginTAD HoxinTAD h |perl -ane 'if($F[2] >0){$sum+=$F[1];};END{print "$sum\n";};'
#13302
#
#
### how many TADS with hox peaks
#cut -f1-3 HoxinTAD |sort -u |wc -l
#3369
#
### now fuse hox in tad and cluster4 in tad
#
### how many TADs have a cluster4 gene ?
#perl fuse.pl cluster_4_up_in_spdh.bed cluster4inTAD HoxinTAD g|cut -f1,1 |sort |uniq -c |wc -l
#49
### how many cluster 4 genes in those TADs ?
#perl fuse.pl cluster_4_up_in_spdh.bed cluster4inTAD HoxinTAD g |perl -ane 'if($F[2] >0){$sum+=$F[2];};END{print "$sum\n";};'
#51
### all those 49 have a hoxd13 peak also in the tad
#
### thus 
#18 genes are not in a TAD at all but may have a hoxd13 peak though, we dont know
#
###### ok, figure 1 
## how many tads with genes
perl fuse.pl $UNG CtrlginTAD HoxinTAD g > tads_with_genes_ctrl
wc -l tads_with_genes_ctrl

## how many tads with genes and hoxd peak
perl fuse.pl $UNG CtrlginTAD HoxinTAD g|perl -ane 'print if($F[1] >0);' > tads_with_genes_and_hoxd13_ctrl
wc -l tads_with_genes_and_hoxd13_ctrl
#2755


## how many tads with genes
perl fuse.pl $C4 cluster4inTAD HoxinTAD g |wc -l
69
## how many tads with genes and hoxd peak
perl fuse.pl $C4 cluster4inTAD HoxinTAD g|perl -ane 'print if($F[1] >0);'|wc -l
67

### using all TADS we get 82 and 80 out of it thus all good!



a=c(80/82,2755/3563)
png("ratio_of_TADS.png")
barplot(a,las=1,yaxt="none",col=c("brown","darkblue"),names.arg=c("Cluster4 genes","Non-cluster4 genes"),ylab="ratio of TADS with gene and hoxd13 peak",ylim=c(0,1))
axis(2,seq(0,1,0.1),paste(seq(0,100,10),"%",sep=""),las=1,ylab="ratio of TADS with gene and hoxd13 peak")
dev.off()

res=c();for(i in 1:10000){od=sample(1:dim(a)[1]);res=c(res,length(which(a[od[1:82],2] > 0)))}

png("shuffling_results.png",1024,1024)
plot(density(res/82),xlim=c(0.5,1),xlab="ratio [Tads with hoxd13/all Tads]",main="sampling 82 tads a 1000 times from all tads\nthe empirical p-value is 0 since none of the 1000 ratios from shuffling was as good as the real one from cluster 4 (red line)")
abline(v=80/82,col="red",lwd=2)
dev.off()

png("shuffling_results_hist.png",1024,1024)
plot(density(res/82),xlim=c(0.5,1),xlab="ratio [Tads with hoxd13/all Tads]",main="sampling 82 tads a 1000 times from all tads\nthe empirical p-value is 0 since none of the 1000 ratios from shuffling was as good as the real one from cluster 4 (red line)")
abline(v=80/82,col="red",lwd=2)
dev.off()


### final final
svglite("shuffling_results_hist.svg",5,10)
hist(res/82,xlim=c(0.5,1),xlab="ratio [Tads with hoxd13/all Tads]",main="sampling 82 tads a 1000 times from all tads\nthe empirical p-value is 0 since none of the 1000 ratios from shuffling\n was as good as the real one from cluster 4 (blue line)",las=1,lwd=2,breaks=10,freq=T,ylim=c(0,500),col="brown")
abline(v=80/82,col="#b7e7ff",lwd=2)
dev.off()

svglite("shuffling_results_hist_density.svg",5,10)
hist(res/82,xlim=c(0.5,1),xlab="ratio [Tads with hoxd13/all Tads]",main="sampling 82 tads a 1000 times from all tads\nthe empirical p-value is 0 since none of the 1000 ratios from shuffling\n was as good as the real one from cluster 4 (blue line)",las=1,lwd=2,breaks=10,freq=F,ylim=c(0,10),col="brown")
abline(v=80/82,col="#b7e7ff",lwd=2)
lines(density(res/82),col="orange")
dev.off()

## final final 3 #### final
svglite("shuffling_results_curve.svg",5,10)
plot(density(res/82),xlim=c(0.5,1),xlab="ratio [Tads with hoxd13/all Tads]",main="sampling 82 tads a 1000 times from all tads\nthe empirical p-value is 0 since none of the 1000 ratios from shuffling\n was as good as the real one from cluster 4 (blue line)",yaxt="n",ylim=c(0,10),col="orange",ylab="Frequency")
abline(v=80/82,col="#b7e7ff",lwd=2)
axis(2,seq(0,10,2),seq(0,500,100),las=1)
dev.off()



## number of hoxd13 peaks in tads
a2=a[which(a[,2] >0),]

## not sure here yet if we do the right thing though

all4 <-read.table("TADS_with_hoxd13_cluster4")
allo <-read.table("TADS_with_hoxd13_all_other")

res=c();for(i in 1:1000){od=sample(1:dim(a2)[1]);res=c(res,mean(a2[od[1:80],2]))}
res2=c();for(i in 1:1000){od=sample(1:dim(a2)[1]);res2=c(res2,a2[od[1:80],2])}

nl <- list(all4[,2],res2)
names(nl)=c("cluster4","1000 times 80\nshuffled TADs")

png("average_hoxd13_per_TAD_and_shuffled.png",1024,1024)
myv(nl,ylims=c(0,40),ttitle="Average Hoxd13 peaks per TAD")
dev.off()


png("average_hoxd13_per_TAD_and_shuffledi_boxplot.png",400,1200)
svglite("average_hoxd13_per_TAD_and_shuffledi_boxplot.svg",4,12)
par(mar=c(8,5,2,1),bty="n",cex=1.2)
boxplot(list(all4[,2],res),las=1,xlab="",xaxt="n",main="Average Hoxd13 peaks\n per TAD",ylab="Hoxd13 peaks per TAD",col=c("cornflowerblue","brown"),bty="n",ylim=c(0,50))
text(x=c(1,2),y=-1.2, srt=45, adj=1, xpd=TRUE,labels=c("TADs with cluster4\n genes (80)","1000 shuffled controls"))


dev.off()


## now figure3 
dist4 <- read.table("closest_hox.csv")
dista <- read.table("closest_hox_allo.csv")

res=c();for(i in 1:10000){od=sample(1:dim(dista)[1]);res=c(res,dista[od[1:84],2])}




png("tss_to_hoxd13_in_TAD_shuffled_version.png",1024,1024)
plot(density(log2(dist4[,2])),ylim=c(0,0.25),col="brown",xlab="log2 distance [nt]",main="cluster4 genes vs Ctrl genes, distance TSS to next hoxd13 peak")
lines(density(log2(res)),col="darkblue")
legend("topleft",c("cluster4","10,000x 84 genes from Ctrl"),lwd=2,col=c("brown","darkblue"))
dev.off()

png("barplot_tss_mean_dist_closest_hox.png",1024,1024)
barplot(c(1785,16384),col=c("brown","darkblue"),names.arg=c("cluster4","10000x 84\nshuffled genes"),ylab="mean distance of TSS to closest hoxd13 peak")
dev.off()




## ok 
a1mb=c(67/69,2497/3265)
png("ratio_of_TADS1mb.png")
barplot(a1mb,las=1,yaxt="none",col=c("brown","darkblue"),names.arg=c("Cluster4 genes","Non-cluster4 genes"),ylab="ratio of TADS with gene and hoxd13 peak",ylim=c(0,1))
axis(2,seq(0,1,0.1),paste(seq(0,100,10),"%",sep=""),las=1,ylab="ratio of TADS with gene and hoxd13 peak")
dev.off()

> res=c();for(i in 1:10000){od=sample(1:dim(a)[1]);res=c(res,length(which(a[od[1:69],2] > 0)))}
> mean(67/69)
[1] 0.9710145
> mean(res/69)
[1] 0.7645855




# figure 2
# how many hoxd13 peaks in tads with cluster4 genes
perl fuse.pl cluster_4_up_in_spdh.bed cluster4inTAD HoxinTAD h |perl -ane 'if($F[2] >0){$sum+=$F[1];};END{print "$sum\n";};'
596

## how many hoxd13 peaks in those tads
perl fuse.pl cluster_4_up_in_spdh.bed cluster4inTAD HoxinTAD h |perl -ane 'if($F[2] >0){ print;}' > TADS_with_hoxd13_cluster4

## how many hoxd13 peaks in those all other Tasds
perl fuse.pl $UNG Ctrlg4inTAD HoxinTAD h |perl -ane 'if($F[2] >0){ print;}' > TADS_with_hoxd13_all_other
mv closest_hox.csv closest_hox_allo.csv

allo <-read.table("TADS_with_hoxd13_all_other")
all4 <-read.table("TADS_with_hoxd13_cluster4")
nl <- list(all4[,2],allo[,2])
names(nl)=c("cluster4 genes","ctrl genes")

png("average_hoxd13_per_TAD.png")
myv(nl,ylims=c(0,40),ttitle="Average Hoxd13 peaks per TAD")

dev.off()


# figure 3
## distance to nearest hoxd13 peak is then what ?
dist4 <- read.table("closest_hox.csv")
dista <- read.table("closest_hox_allo.csv")


png("tss_to_hoxd13_in_TAD.png")
plot(density(log2(dist4[,2])),ylim=c(0,0.25),col="brown",xlab="log2 distance [nt]",main="cluster4 genes vs Ctrl genes, distance TSS to next hoxd13 peak")
lines(density(log2(dista[,2])),col="darkblue")
legend("topleft",c("cluster4","Ctrl"),lwd=2,col=c("brown","darkblue"))
dev.off()





############# better violin plots
myv <- function(x,ttitle="violinPlot",ylims=c(-4,4)){
    require(ggplot2)
    categories=length(x)
    r.vals=unlist(x)
    indicator = rep(categories,length(r.vals))
    start=1
    for(k in 1:categories){
        indicator[start:(start+length(x[[k]])-1)]=k
        start=start+length(x[[k]])

    }
    tt=factor(indicator,ordered=T,labels=names(x))
    fmat2=data.frame(values=r.vals,type=as.numeric(tt),supp=tt)
    p<-ggplot(fmat2,aes(x=supp, y=values, fill=supp))+geom_violin(trim=T)+
        theme_classic(base_size=20)+
        labs(x="",y="average hoxd13 peaks per TAD",fill="type",title=ttitle)+
        scale_fill_brewer(palette="cornflowerblue")+
        ylim(ylims[1],ylims[2])

    nlabels <- unlist(lapply(x,length))
    nnlabels<-paste(labels(nlabels),paste("(",nn(as.numeric(nlabels)),")",sep=""),sep="\n")

    p<-p+ scale_x_discrete(labels= nnlabels)
    return(p)
}




## TSS distances
hist(log2(dist4[,2]),freq=F,breaks=10,col="cornflowerblue",las=1,xlab="Peak distance from TSS [nt]",main="Distribution of distances from TSS to closest Hoxd13 peak",xaxt="n")
axis(1,4:18,prettyNum(2^(4:18)))
lines(density(log2(dist4[,2])),col="red")



average  ....
boxplots


### #### mmm 
svglite("tss_distance.svg",20,7)
pdf("tss_distance.pdf",20,5)

par(mfrow=c(1,2))

hist(log2(dist4[,2]),freq=F,breaks=10,col="cornflowerblue",las=1,xlab="Peak distance from TSS [nt]",main="Distribution of distances from TSS to closest Hoxd13 peaks\nCluster 4 dysregulated genes",xaxt="n",ylim=c(0,0.2))
axis(1,4:18,prettyNum(2^(4:18),big.mark=","))
lines(density(log2(dist4[,2])),col="red")
abline(mean=)


hist(log2(res),freq=F,col="brown",las=1,xlab="Peak distance from TSS [nt]",main="Distribution of distances from TSS to closest Hoxd13 peak\n1000 times random selected genes",xaxt="n",breaks=20)
axis(1,seq(0,23,2),prettyNum(2^seq(0,23,2),big.mark=","))
lines(density(log2(dist4[,2])),col="orange")
dev.off()


dist4 <- read.table("closest_hox.csv")
dista <- read.table("closest_hox_allo.csv")

resX=c();for(i in 1:1000){od=sample(1:dim(dista)[1]);res=c(res,dista[od[1:84],2])}
wt(resX)


resX=read.table("resX")


## final final
for(fin in c(T,F)){
    if(fin == T){
        svglite("tss_distance_freq.svg",20,7)
        ylims1=c(0,14)
        ylims2=c(0,2e4)
    }else{
        svglite("tss_distance_density.svg",20,7)
        ylims1=c(0,0.25)
        ylims2=c(0,0.25)
    }
#pdf("tss_distance_freq.pdf",20,7)
par(mfrow=c(1,2))
xran <- c(0,24)
hist(log2(dist4[,2]),freq=fin,breaks=10,col="cornflowerblue",las=1,xlab="Peak distance from TSS [nt]",main="Distribution of distances from TSS to closest Hoxd13 peaks\nCl
uster 4 dysregulated genes",xaxt="n",ylim=ylims1,xlim=xran)
ran=seq(0,24,2)
ggo<-c(2^ran)
ggg<-round(ggo/1024)
xlabx=c(ggo[which(ggg==0)], paste(ggg[which(ggg < 1000 & ggg > 0)],"k",sep=""),paste(round(ggg[which(ggg > 1000)])/1024,"M",sep=""))
axis(1,ran,xlabx)
#if(length(which(ggg > 1000)) > 0){
#   xlabx=c(xlabx,paste(round(ggg[which(ggg > 1000)])/1024,"M",sep=""))
#}
#lines(density(log2(dist4[,2])),col="orange")
abline(v=log2(mean(dist4[,2])),col="#b7e7ff",lwd=2)
hist(log2(resX),freq=fin,col="brown",las=1,xlab="Peak distance from TSS [nt]",main="Distribution of distances from TSS to closest Hoxd13 peak\n1000 times random selected
genes",xaxt="n",breaks=20,ylim=ylims2,xlim=xran)
axis(1,ran,xlabx)
d2<-density(log2(resX))
to2<-which(d2$x < 1)
#lines(d2$x[-to2],d2$y[-to2],col="orange")
abline(v=log2(mean(resX)),col="#b7e7ff",lwd=2)
dev.off()
}










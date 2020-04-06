#!/usr/local/bin/Rscript

# generate_figures_S7A-E.R -  R script
# Copyright (C) 2018  Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

## here are lying general R scripts
DPATH="../../Fdata/"
SPATH="../../src/"

library(dendextend)
library(svglite)
library(gplots)
source(paste(SPATH,"colors.R",sep=""))
source(paste(SPATH,"functions.R",sep=""))
source("../../Figure7/Figure7A/labbook_analysis_functions.R")
source("../../Figure7/Figure7A/circplot_fun_re.R_new")

## adapt so it an be easily used by all 
pos_set=read.table(paste(DPATH,"Figure7/lambert/ENSG_TRUE_tf_lambert_zf_split_AA_repeat_lengths_8plus.csv",sep=""))
pos_set10=read.table(paste(DPATH,"Figure7/lambert/ENSG_TRUE_tf_lambert_zf_split_AA_repeat_lengths_10plus.csv",sep=""))

DBD=paste(DPATH,"Figure7/lambert/2018_AddAlignments/",sep="")

## all sets now length filtered
tf.list  <- readAAm(sel=1:32,insets=tfsets)
dbd.list <- readAAm(sel=1:4,insets=dbdsets,P=DBD,addf='')

pmat.tf2_10 <- process.mat(tf.list,"sAA",thresAA=0.10)
pmat.dbd_10 <- process.mat(dbd.list,"sAA",thresAA=0.10) ## we have many less TF left with single AA enrichment for DBDs

tlax<-circplot(pmat=pmat.tf2_10,tfl=tf.list,pos_set=pos_set10,nfam=8,nclust=7,dims.use=10)
dev.off()

## we need the new gene order here in a vector so we can use it for figure S7E
## the plot is the one in figure S7A but we don't need it here
new_gene_order<-replot_circle(tlax,plotdend=F,cexl=1.4,colsFam=j.colors10,colsClust=colors7[c(1,4,3,5,6,2,7)])
dev.off()

### Figure S7A
dmat=pmat.tf2_10[[1]]
cmat=pmat.dbd_10[[1]]
bmat=rbind(dmat,cmat)

b_dall <-getdend(bmat,dims.use=10,nclust=7,md="ward.D2")
df1 <- hkm(b_dall,both=T)

gene_ids_cluster=as.numeric(b_dall[[4]]$cluster)

## we make a circle plot now out of it!!! They are much simpler than the IDR only plot so far!!!! hopefully 
df1_labels <- as.numeric(labels(df1))             ## these are the labels in linear order from left to right -> we will plot it like this
df1_gene_clusters <- gene_ids_cluster[df1_labels]
nsize=length(df1_labels) ## just the length

dbds.i<-which(df1_labels > nrow(dmat))
idrs.i<-which(df1_labels <= nrow(dmat))


idrs.col=rep(colors7[3],7)
dbds.col=rep(colors7[6],7)

mytcol=rep("white",nsize)

mytcol[dbds.i]=dbds.col[df1_gene_clusters[dbds.i]]
mytcol[idrs.i]=idrs.col[df1_gene_clusters[idrs.i]]

doPng=1
if(doPng){
	png("FigureS7A.png",1000,1000)
}else{
	svglite(file = "FigureS7A.svg", width = 12,height=12)
}

circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.00),gap.degree=0)
circos.initialize("foo", xlim = c(0, nsize))
mtrack(mycol=mytcol,th=0.06)


legend("topright",c("Intrinsically disorde-\nred regions (IDRs)"," ","DNA binding do-\nmaines (DBDs)"),col=c(idrs.col[1],"white",dbds.col[1]),pch=15,border="grey")

dend=df1
max_height = attr(dend, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
circos.dendrogram(dend, max_height = max_height)}, track.height = 0.4, bg.border = NA)
text(0,0,"IDRs + DBDs\n(1446)  (365)")
dev.off()


### Figure S7B
cmat=pmat.dbd_10[[1]]
dbd_dall <-getdend(cmat,dims.use=10,nclust=4,md="ward.D2",seed=8818)

## check again
df1 <- hkm(dbd_dall,nclust=4)
df1_labels <- as.numeric(labels(df1))             ## these are the labels in linear order from left to right -> we will plot it like this
gene_ids_cluster <- as.numeric(dbd_dall[[4]]$cluster)
df1_gene_clusters <- gene_ids_cluster[df1_labels]
nsize=length(df1_labels) ## just the length

## orig index
famind <- dbd.list$groupindex[pmat.dbd_10$indices]
## reordered
famind.r <- famind[df1_labels]

dbd.col.only=colors7[c(1,6,3,4)] ## corresponds now to TF families 1,2,3,4
dbd.col.clust=dbd.col.only[c(4,1,3,2)]



for(tt in c("png","svg")){
if(tt == "png"){
png("FigureS7B.png",1000,1000)
}else{
svglite(file = "FigureS7B.svg", width = 12,height=12)
}

circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.06),gap.degree=0)
circos.initialize("foo", xlim = c(0, nsize))

mtrack(mycol=dbd.col.only[famind.r],th=0.04)
legend("topright",legend=dbd.list$sets,col=dbd.col.only[1:4],pch=16)

## the clusters
mtrack(mycol=dbd.col.clust[labels_col(df1)],th=0.04)
#legend(-1,1,legend=1:4,col=dbd.col.clust[1:4],pch=16)

dend=df1
max_height = attr(dend, "height")
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.00),gap.degree=0)
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
circos.dendrogram(dend, max_height = max_height)}, track.height = 0.5, bg.border = NA)
text(0,0,"DBDs\n(365)")

dev.off()
}

#################################################
### Figure S7C
#png("FigureS7C.png",800,800)
svg("FigureS7C.svg",8,8)
myscree(tlax$idrs_dall$pca,pca=T,cs=T,varline=80)
dev.off()

#################################################
## Figure S7D

#bicplot.fit(tlax$idrs_dall$pca,s=999)

tmp <- bicplot(tlax$idrs_dall$pca,seed=999,biconly=T,plotme=F)
a=matrix(unlist(tmp)[seq(2,63,3)])

svg("FigureS7D.svg",12,12)
plot(1:21,a[,1],"l",col="cornflowerblue",xaxt="n",las=2,ylab="Bayesian information criterion",xlab="k-means clusters" ,bty="n",main="")
axis(1,1:22,2:23)
abline(v=6,col="lightgrey")
abline(h=a[6,1],col="lightgrey")
dev.off()

## Figure S7E
feature.mat=matrix(0,nrow=7,ncol=25)
colnames(feature.mat)=colnames(pmat.tf2_10$matrix)
rownames(feature.mat)=paste("C",1:7)
fm=feature.mat
aalist <- read.table("aa.csv",header=F,sep="\t")

coln<- colnames(feature.mat)
get_n <-rep(0,length(coln))
for(i in 1:length(coln)){
    res <- which(aalist[,3] %in%coln[i])
    if(length(res) >0){
        get_n[i] = res
    }
}
get_n
coln2 <- coln
torep <- which(get_n > 0)
coln2[torep] <- paste(addbraces(aalist[,3]),as.vector(aalist[,1]))

colnames(feature.mat)<-coln2
fm=feature.mat

## KS-testing
mm <- fm
mm.l <- fm
cellnotes=mm
obj <- tlax$idrs_dall
obj[[6]] <- pmat.tf2_10$matrix
names(obj)=c(names(tlax$idrs_dall),"matrix")
for(param in 1:length(coln)){
#for(param in 1){
    ## reference cluster
    for(ci in 1:7){
        ## get the new gene ids here we have from the replot function
        cin=which(new_gene_order==ci)
        ## density
        da=density(obj[[6]][cin,param])
        cin.ymax<-which(da$y==max(da$y))
        cellnotes[ci,param]=round(da$x[cin.ymax],2)
        #cellnotes[ci,param]=da$x[cin.ymax]

        ## get all other clusters
        ni <- c(1:7)[-ci]

        ## pool rest and ask for ci and we take the max of the density peak here !!!
        ks <- ks.test(obj[[6]][-cin,param],obj[[6]][cin,param],alternative="greater")
        ## depletion
        ks2 <- ks.test(obj[[6]][-cin,param],obj[[6]][cin,param],alternative="less")

        if(ks$p.value < 0.01/(7*25)){
            mm[ci,param] <- -log10(ks$p.value)
        }
        if(ks2$p.value < 0.01/(7*25)){
            if(mm[ci,param] == 0){
                log10(ks2$p.value)
            }else{
                print(paste(param,"in cluster",ci,"is enriched",mm[ci,param]))
            }
            mm[ci,param]<- log10(ks2$p.value)
        }
        ## get peak density of param
    }
}


## set upper and lower bound for plotting the heatmap
mm2=mm
mm2[mm2>50]=51
mm2[mm2< (-50)]= (-51)

## order by cluster 6, the 
no=order(mm2[6,],decreasing=T)
mm3=mm2[,no]
ncsa=order(colSums(mm2),decreasing=T)

col.here=colorRampPalette(c("blue","white","red"))(n=101)

svg("FigureS7E.svg",20,8)
par(oma=c(6,5,1,1))
heatmap.2(mm3[,ncsa],dendrogram="none",trace="none",
		  Rowv=F,Colv=T,
          col=col.here,las=2,
          srtCol=45,key=T,RowSideColors=colors7,
          key.xlab="log10 adj. p-values\n(ks.test)",
          density.info="histogram",keysize=1.1,cexCol=1.6,offsetRow=-121,cexRow=2,
          ## no cell notes
          #cellnote=signif(mm3[,ncsa],4),notecol="black",notecex=1.5,
          key.title="depleted              enriched")
dev.off()

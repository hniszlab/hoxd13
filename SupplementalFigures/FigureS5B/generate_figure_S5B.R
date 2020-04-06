#!/usr/bin/Rscript

# generate_figure_S5B.R -  R script
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

DPATH="../../Fdata/Figure5B/"
SPATH="../../src/"
source(paste(SPATH,"colors.R",sep=""))

library(Matrix)
library(Rtsne)

## package for fast SVD and thus PCA 
library(irlba)

## load package for additional functions to visualize
library(Seurat,lib.loc=paste(DPATH,"/myseurat2.3.4",sep=""))

source("../../Figure5/ScRNA-Seq_data_processing/process_seurat.R")

#adjust file names depending on cell ranger output or just use your own ones
countmatrix="matrix.mtx"
rownames="features_u.tsv"
colnames="barcodes.tsv"

load10x <- function(p="./",rown=rownames,coln=colnames,matn=countmatrix){
	mat <- readMM(paste(p,"/",matn,sep=""))
	rows <- read.table(paste(p,"/",rown,sep=""),sep="\t")
	cols <- scan(paste(p,"/",coln,sep=""),what="char")
	colnames(mat)=cols
	rownames(mat)=rows[,2]
	return(mat)
}


cwd=getwd()

mysample ="wt"
setwd(paste(DPATH,mysample,sep=""))
ctrl <- load10x("./")
setwd(cwd)

mysample="spdh"
setwd(paste(DPATH,mysample,sep=""))
spdh <- load10x("./")
setwd(cwd)


## which cells are hoxd13 positive
mysample="wt"
p=paste(DPATH,mysample,sep="")
rows <- read.table(paste(p,"/",rownames,sep=""),sep="\t")



## get number of cells in both conditions
cells_spdh <- dim(spdh)[2]
cells_ctrl<- dim (ctrl)[2]

## just getting a limiting number of cells here
bord <- 5000
if(cells_spdh < 1e4){
	bord<- 1e3*ceiling(cells_spdh/1e3)
}else{
	bord<- 1e4*ceiling(cells_spdh/1e4)
}


## set colnames appropriate
colnames(spdh)=paste("cell",1:cells_spdh,sep="_")
##                           5001     5000 
colnames(ctrl)=paste("cell",(bord+1):(bord+cells_ctrl),sep="_")

## fuse both matrices
fullmat<-cbind(spdh,ctrl)

## these are the gene symbols only 
rw=rownames(fullmat)
hoxd13=which(rw=="Hoxd13")

## which are active and which not
hoxd13_active=which(fullmat[hoxd13,] > 0)
hoxd13_inactive=which(fullmat[hoxd13,] == 0)


## make seurat object now
pbmc <- CreateSeuratObject(raw.data = fullmat[,hoxd13_active], min.cells = 3, min.genes = 200, project = "limb")

## get genes on chromosome M
mito=scan(paste(DPATH,"gencode.vM19_mitochondrial_genes.id",sep=""),what="char")
pbmc_no_regress <- process(pbmc)

## I can first plot the cell fate map 
t1t=pbmc_no_regress@dr$tsne@cell.embeddings
p=pbmc_no_regress
rnnt <- as.numeric(unlist(strsplit(rownames(t1t),"_"))[seq(2,2*dim(t1t)[1],2)])

spdh.it <- which(rnnt < bord)
ctrl.it <- which(rnnt >= bord)


mycolors=c("#a65628","#c49a6c","#add8e6","#00008b","#377eb8","#ff7f00","#fbb040","#f9ed32","#4daf4a","#f881bf","#e41a1c")[c(1,2,3,9,5,6,7,8,4,10,11)]

#barplot(pvars,main="Eigenwerte aka squared sd = var per PC")
png("cluster_hoxd13_positive_fate_spdh_wt_no_regress_neu.png",1800,600)

#library(svglite)
#svglite("cluster_hoxd13_positive_fate_spdh_wt_no_regress_neu.svg",18,6)
i=1
j=2
par(mfrow=c(1,3))

plot(t1t[,1],t1t[,2],xlab="t-SNE 1",ylab="t-SNE 2",main="",col="white",las=1)
for(k in 0:8){
	rk=k+1
    points(t1t[p@ident==k,1],t1t[p@ident==k,2],col=mycolors[rk],pch=16)
}
legend("topright",legend=1:9,col=mycolors[1:9],pch=16,bty="n",cex=1.2,xjust=1)

plot(t1t[,i],t1t[,j],xlab="t-SNE 1",ylab="t-SNE 2",main="",col="white")
points(t1t[ctrl.it,i],t1t[ctrl.it,j],col="blue",pch=16)
legend("topleft",c(paste("ctrl",length(ctrl.it))),col=c("blue"),pch=15)

plot(t1t[,i],t1t[,j],xlab="t-SNE 1",ylab="t-SNE 2",main="",col="white")
points(t1t[spdh.it,i],t1t[spdh.it,j],col="brown",pch=16)
legend("topleft",c(paste("spdh",length(spdh.it))),col=c("brown"),pch=15)
dev.off()

## table for S5C
## cluster8 is counted from 0, so it is cluster9 in the figure
markers_no_regress <- FindAllMarkers(object = pbmc_no_regress, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
cluster8 <- markers_no_regress[which(markers_no_regress$cluster == 8),]
write.table(cluster8[,7],file="cluster8_genes.tsv",col.names=F,row.names=F,quote=F)





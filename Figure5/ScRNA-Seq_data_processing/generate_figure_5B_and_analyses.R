#!/usr/local/bin/Rscript

# generate_figure_5B_and_analyses.R -  R script
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

library(svglite)
library(Matrix)
library(Rtsne)
library(irlba)

## the script was created with Seurat version 2.3.4.
library(Seurat,lib.loc=paste(DPATH,"myseurat2.3.4",sep=""))

## get genes on chromosome M
mito=scan(paste(DPATH,"/gencode.vM19_mitochondrial_genes.id",sep=""),what="char")

source("./process_seurat.R")

#adjust file names depending on cell ranger output or just use your own ones
countmatrix="matrix.mtx"
rownames="features_u.tsv"
colnames="barcodes.tsv"

## loading function
rows <-c()
load10x <- function(p="./",rown=rownames,coln=colnames,matn=countmatrix){
	mat <- readMM(paste(p,"/",matn,sep=""))
	rows <<- read.table(paste(p,"/",rown,sep=""),sep="\t")
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


## just getting a number here as an index discerning both populations
bord <- 5000
if(cells_spdh < 1e4){
	bord<- 1e3*ceiling(cells_spdh/1e3)
}else{
	bord<- 1e4*ceiling(cells_spdh/1e4)
}

## keep the original cell names aka barcode used since we need it for the browser tracks

## uncomment if you want the cell number instead
## set colnames appropriately
#colnames(spdh)=paste("spdh",1:cells_spdh,sep="_")
##                           5001     5000 
#colnames(ctrl)=paste("ctrl",(bord+1):(bord+cells_ctrl),sep="_")

## these are the gene symbols only 
rw=rownames(ctrl)
hoxd13=which(rw=="Hoxd13")
## which are active and which not
hoxd13_active=which(ctrl[hoxd13,] > 0)
hoxd13_inactive=which(ctrl[hoxd13,] == 0)


## make seurat object now for all control runs
## we take all cells not just the positive ones 
pbmc.ctrl <- CreateSeuratObject(raw.data = ctrl,min.cells = 3, min.genes = 200, project = "limb")
#pbmc.act <- CreateSeuratObject(raw.data = ctrl[,hoxd13_active], min.cells = 3, min.genes = 200, project = "limb")
#pbmc.ina <- CreateSeuratObject(raw.data = ctrl[,hoxd13_inactive], min.cells = 3, min.genes = 200, project = "limb")

## do the processing now. The process function is in script ./process_seurat.R
pbmc.ctrl <- process(pbmc.ctrl)
#pbmc.act <- process(pbmc.act)
#pbmc.ina <- process(pbmc.ina)


hoxd13.s=which(rownames(spdh)=="Hoxd13")
#hoxd13_active.s=which(spdh[hoxd13.s,] > 0)
#hoxd13_inactive.s=which(spdh[hoxd13.s,] == 0)


## make all spdh runs
pbmc.spdh <- CreateSeuratObject(raw.data = spdh ,min.cells = 3, min.genes = 200, project = "limb")
#pbmc.spdh.act <- CreateSeuratObject(raw.data = spdh[,hoxd13_active.s] ,min.cells = 3, min.genes = 200, project = "limb")
#pbmc.spdh.ina <- CreateSeuratObject(raw.data = spdh[,hoxd13_inactive.s] ,min.cells = 3, min.genes = 200, project = "limb")


pbmc.spdh <- process(pbmc.spdh)
#pbmc.spdh.act <- process(pbmc.spdh.act)
#pbmc.spdh.ina <- process(pbmc.spdh.ina)


## put all in one list
#allc <-list(pbmc.ctrl,pbmc.act,pbmc.ina,pbmc.spdh,pbmc.spdh.act ,pbmc.spdh.ina)
allc <-list(pbmc.ctrl,pbmc.spdh)
names(allc) <- c("ctrl","spdh")
#names(allc) <- c("ctrl","ctrl_hox_pos","ctrl_hox_neg","spdh","spdh_hox_pos","spdh_hox_neg")

## raw plot from figure 5B
#TSNEPlot(object = allc[[1]],plot.title=names(allc)[1])

## plot in figure 5B
p <- pbmc.ctrl
current.cluster.ids=0:10
new.cluster.ids=1:11
p@ident<- plyr::mapvalues(p@ident, from = current.cluster.ids, to = new.cluster.ids)
t1t=p@dr$tsne@cell.embeddings
mycolors=c("#a65628","#c49a6c","#add8e6","#00008b","#377eb8","#ff7f00","#fbb040","#f9ed32","#4daf4a","#f881bf","#e41a1c")
mydims=c(1,2) ## not sure why it said 3 before. we only have 2d in tsne

#############################################################################################
#
#
# make Figure 5B now
#
#
#############################################################################################

## plot it 
#for(i in 1:6){
#	x11()
#	TSNEPlot(object = allc[[i]],plot.title=names(allc)[i])
#}

library(svglite)
svglite("Figure5B.svg",15,15)
par(cex=1.2)
plot(t1t[,1],t1t[,2],xlab="t-SNE 1",ylab="t-SNE 2",main="",col="white",xlim=c(-45,50),las=1)
for(k in 1:11){
    points(t1t[p@ident==k,mydims[1]],t1t[p@ident==k,mydims[2]],col=mycolors[k],pch=16)
}
legend("topright",legend=1:11,col=mycolors,pch=16,bty="n",cex=1.2,xjust=1)
dev.off()

###################################################
#
#
# Marker gene identification
#
#
###################################################


## find marker genes now being differentially expressed between the clusters
markers = list()
#for(i in 1:6){
for(i in 1){
	markers[[i]] <- FindAllMarkers(object = allc[[i]], only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) 
}


##################################################
#
#
# Average expression of genes in wt clusters 
#
#
##################################################

## get average expression for all genes in all clusters
## this is exp( rpm(umi)) -1 
## this actually reverts the Logartihm in the preprocessing but not the Normalization.

## this will have 11 columns, one for each cluster right now counting from 0 to 11
avex=list()
#for(i in 1:6){
for(i in 1){
	avex[[i]] <-AverageExpression(object = allc[[i]])
}

i=1
sobj=pbmc.spdh@data

## check first which genes are present in both conditions 
gin <- which(rownames(avex[[i]]) %in% rownames(sobj))
## same 11 clusters counting from 0 to 10
## these are now the ctrl average values for the set of genes we also find in the spdh object
leftavex <- avex[[i]][gin,]

## now the reverse
gin2 <- which(rownames(sobj)%in%rownames(leftavex))
spdh_to_assign <- sobj[gin2,]

## both objects have the same length
length(gin)
length(gin2)


## an now we reorder sobj so we get the same order as in leftavex
nord=c()
rna <- rownames(leftavex)
rns <- rownames(spdh_to_assign)

for(i in 1:length(gin2)){
	nord=c(nord,which( rns %in% rna[i]))
}

## reorder if necessary
spdh_avex=spdh_to_assign[nord,]

## put in same dimensional scale as average cluster values are per gene, so we need to do expm1 on the spdh@data object
spdh_avex_expm1=expm1(spdh_avex)


##################################################
#
#
# Assignment of SPDH cells to wt clusters
#
#
##################################################

## euclidian distance
dfun <- function(x,y){
	sqrt(sum((x-y)^2))
}

## assign cell to nearest cluste
spdh_avex.cl=c()
for(i in 1:dim(spdh_avex_expm1)[2]){
	spdh_avex.cl=c(spdh_avex.cl,order(apply(leftavex,2,dfun,y=spdh_avex_expm1[,i]))[1])
}

### wt cluster cell numbers
table(pbmc.ctrl@ident)
#  1   2   3   4   5   6   7   8   9  10  11 
#722 614 593 564 564 540 410 301  70  59  27 

## spdh_avex.cl -> goes from 1:11 mapping 0:10 in leftavex, since order and sort are increasing the first element is the smallest 
## which is the one with smalles distance
table(spdh_avex.cl)
#  1   2   3   4   5   6   7   8   9  10  11 
#765 415 433 270 637 867 269 348  57  39  47 

exprs_p<-list()
for(i in 0:10){
	j=i+1
	names_markers <- markers[[1]][which(markers[[1]]$cluster==i),7]
	indexm <- which(rownames(avex[[1]])%in%names_markers)

	exprs_p[[j]] <- as.numeric(avex[[1]][indexm,j])
	exprs_p[[j]] <-cbind(exprs_p[[j]],rowMeans(avex[[1]][indexm,setdiff(1:11,j)]))
	exprs_p[[j]] <-cbind(exprs_p[[j]],log(exprs_p[[j]][,1]/exprs_p[[j]][,2]))

	rownames(exprs_p[[j]]) <- rownames(avex[[1]][which(rownames(avex[[1]])%in%names_markers),])
	colnames(exprs_p[[j]])=c("exprs_marker_cluster","avg_exprs_others","logFC")
}


## these are the actual mean expression values for the marker genes in the 11 different clusters 
## Printing
if(!file.exists("marker_genes")){
    dir.create(file.path("marker_genes"))
}

for (i in 0:10){
	j=i+1
	write.table(exprs_p[[j]],file=paste("marker_genes/mean_exprs_ctrl_clusterN_",j,".tsv",sep=""),quote=F)
}


###############################################################
#
#
# Differential expression analysis
#
#
###############################################################

## combine both ctrl and spdh seurat objects
pbmc.combined <- MergeSeurat(object1 = pbmc.ctrl, object2 = pbmc.spdh, add.cell.id1 = "ctrl",add.cell.id2 = "spdh",project="WTSPDH")
## process them
pbmc.combined <-process(pbmc.combined)

## spdh cells to assigned clusters from ctrl 1-11 get cluster ids 21-31 so we can do a diffex analysis  
nnames=as.factor(c(as.numeric(unlist(pbmc.ctrl@ident)),spdh_avex.cl+20))

## name the nnames cells
#names(nnames)=c(paste("ctrl_",colnames(pbmc.ctrl@data),sep=""),paste("spdh_",colnames(sspd@data),sep=""))
names(nnames)=c(paste("ctrl_",colnames(pbmc.ctrl@data),sep=""),paste("spdh_",colnames(pbmc.spdh@data),sep=""))

## assign new names to them 
pbmc.combined@ident=nnames

## and now we compare the cells in the ctrl cluster i=1-11 to the spdh cells in the corresponding spdh cluster 21-31
## the FindMarkers function makes diffex testing based on Wilcoxon rank sum test
allm <-list()
for(i in 1:11){
    allm[[i]] <- FindMarkers(pbmc.combined, ident.1 = i, ident.2 = (i+20))
}


## make dir if not existent
if(!file.exists("diffex_genes")){
	dir.create(file.path("diffex_genes"))
}

# the output values used for the test are scaled are thus different from what is gotten from the Average expression function.
### write out data
for (i in 1:11){
    j=i-1
    write.table(allm[[i]], file=paste("diffex_genes/diffex_genes_cluster_",i,".tsv",sep=""),quote=F,row.names=T,col.names=c("p_val","avg_log2FC","wt","spdh","p_val_adj"))
}


####################################
#
#
# Figure5F/S5D need the cell barcodes
#
#
###################################

if(!file.exists("barcodes")){
    dir.create(file.path("barcodes"))
}


## write out barcodes of each cell per cluster for ctrl cells
for(i in 0:10){
    j=i+1
    ttt<-names(pbmc.ctrl@ident)[which(pbmc.ctrl@ident == i)]
    write.table(ttt, file=paste("./barcodes/barcodes_ctrl_cells_cluster_",j,".tsv",sep=""),quote=F,row.names=F,col.names=F)
}

## write out barcodes of each cell per clusters for spdh assigned
for(i in 0:10){
    j=i+1
    ttt<- colnames(sobj)[which(spdh_avex.cl == j)]
    write.table(ttt, file=paste("./barcodes/barcodes_spdh_cells_cluster_",j,".tsv",sep=""),quote=F,row.names=F,col.names=F)
}

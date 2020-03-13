#!/usr/bin/Rscript
library(Matrix)
library(Rtsne)

## package for fast SVD and thus PCA 
library(irlba)

## load package for additional functions to visualize
library(Seurat,lib.loc="../../Figure5/ScRNA-Seq_data_processing/seurat2.3.4")

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
setwd(paste("/project/hnisz_lab_analyses/scRNAseq_limb/",mysample,"/filtered_feature_bc_matrix",sep=""))
ctrl <- load10x("./")

mysample="spdh"
setwd(paste("/project/hnisz_lab_analyses/scRNAseq_limb/",mysample,"/filtered_feature_bc_matrix",sep=""))
spdh <- load10x("./")

## going back to orig dir
setwd(cwd)

## which cells are hoxd13 positive
mysample="wt"
p=paste("/project/hnisz_lab_analyses/scRNAseq_limb/",mysample,"/filtered_feature_bc_matrix",sep="")
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
mito=scan("/project/hnisz_lab_storage/annotation/gencode.vM19_mitochondrial_genes.id",what="char")
if(0){
rows.f<-rownames(pbmc@raw.data)

## get ens_gene_ids that are in the mitochondrial gene fraction
rows.mito.i <- which(rows[,1]%in%mito)

## table with ensgene and gene symbol only
rows.mito<-rows[rows.mito.i,1:2]

## indices in rows all that are still in after filtering
mito.genes <- which(rows.f%in%rows.mito[,2])

percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
## take out cells with less than 200 genes or more then 7000 and more than 5% of mitochondrial reads
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(7000, 0.05))
## normalize now
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
## variable genes
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

###############
###############
pbmc_no_regress <- pbmc
}
pbmc_no_regress <- process(pbmc)

## I can first plot the cell fate map 
t1t=pbmc_no_regress@dr$tsne@cell.embeddings
p=pbmc_no_regress
rnnt <- as.numeric(unlist(strsplit(rownames(t1t),"_"))[seq(2,2*dim(t1t)[1],2)])

spdh.it <- which(rnnt < bord)
ctrl.it <- which(rnnt >= bord)

## scree plot
#pvars=((pbmc@dr$tsne@sdev^2/sum(pbmc@dr$pca@sdev^2))*100)

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


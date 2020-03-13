#!/usr/bin/Rscript

# generate_figure_GWAS.R -  R script
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
source(paste(SPATH,"colors.R",sep=""))

library(svglite)
library(gplots)

outmat=read.table(paste(DPATH,"Figure7/outmatrix_v4.csv")
colsin=foldc.colors9[5:9]
breaksin=c(0,1.30103,3,5,10,50)
csizes=table(outmat[,3])

source("./gwas_to_gene_all.R")
gwas.names=names(gwas) ## all the names of the different phenotypes

GWenrich=list() ## make an output list
GWenrich2=list() ## make an output list
signif=c()
GWenrich2_len=c()
matrixAAin=matrix(0,nrow=7,ncol=length(gwas))

for(p.ind in 1:length(gwas)){
	phenot.genes = gwas[[p.ind]] ##  family genes
	len.pg = length(phenot.genes)
	len2.pg= length(which(phenot.genes %in% outmat[,1]))
	GWenrich2_len[p.ind] = len2.pg
	
	res.activ=c() ## result vector
	res.activ2=c() ## result vector
	for(i in 1:7){ ## foreach cluster
## get the cluster ids
		ci=which(outmat[,3] == i)
		aic<-which(phenot.genes %in% outmat[ci,1]) ## which pheno genes in cluster X

## pheno genes not not in cluster , ## but actually we need to know how many we have really in the dataset !!!! -> adjust len.pg for it
		anic <- len.pg-length(aic)
		anic2 <- len2.pg-length(aic)

## index in csizes that is our cluster
		clust.i<-which(names(csizes) == i)

## what is left in cluster that is not a Activ/Repress
		caleft=csizes[clust.i]-length(aic)

## what is globally left
#		atotal=nrow(outmat)-length(aic)-anic-caleft
		atotal2=nrow(outmat)-length(aic)-anic2-caleft

#		mmat=matrix(c(length(aic),anic,caleft,atotal),ncol=2)
		mmat2=matrix(c(length(aic),anic2,caleft,atotal2),ncol=2)
		matrixAAin[i,p.ind]=length(aic)
#		ft <- fisher.test(mmat,alternative="greater")         ## check precisely. maybe we can use it for enrichment "g" and depletion "l"
		ft2 <- fisher.test(mmat2,alternative="greater")         ## check precisely. maybe we can use it for enrichment "g" and depletion "l"
#		res.activ=c(res.activ,ft[[1]])
		res.activ2=c(res.activ2,ft2[[1]])
	}
#	PHenrich[[p.ind]]=res.activ
	GWenrich2[[p.ind]]=res.activ2
	if(length(which(res.activ2 < 0.05)) > 0){
		signif=c(signif,p.ind)
	}

}


## take only those categories which have more then 9 members in our TF set
#tolook_at <- signif[which(GWenrich2_len[signif]>=15)]
resmat.tmp <-matrix(unlist(GWenrich2),byrow=T,ncol=7)   ## in this case it is true since each list entry in PHenrich is a category which has 7 pvalues
rownames(resmat.tmp)=paste(gwas.names,addbraces(GWenrich2_len))

## transpose matrix
## the term was too long so I use a shrinked version!!! for it
rownames(resmat.tmp)[1832]="Chronic_inflammatory_diseases"
resmat <- t(resmat.tmp)


## filter here a bit
takeout<-which(GWenrich2_len<20)                         ## filter on min set size 
takeout<-c(takeout,which(GWenrich2_len > 1000))             ## filter on max set size

takeout<-c(takeout,which(apply(resmat,2,min) >= 0.05)) ## filter on q/p -value
takeout<-c(takeout,which(apply(matrixAAin,2,max)< 15)) ## filter on abs number of enriched genes per set in cluster  
## order by max value per column
takeout <- unique(takeout)
no <- order(apply(matrixAAin[,-takeout],2,max),decreasing=T)
####

svglite("Figure7D_GWAS.svg",12,10)

#par(oma=c(50,1,1,5))
par(oma=c(30,15,1,5))
ks=1
of1=0
of2=0
heatmap.2(-log10(resmat[,-takeout][,no]),
    dendrogram="none",Rowv=F,Colv=F,
    tracecol=NA,
#    RowSideColors=colors7,
#   xlab="Enrichment for members\ncontaining the indicated\nAA repeat",
#    offsetRow=-of1,offsetCol=-of2,
#    sepcolor="black",sepwidth=c(0.005,0.005),colsep=0:ncol(resmat),rowsep=0:nrow(resmat),
    labRow=paste("IDR cluster",1:7),
	breaks=breaksin,
	col=colsin,
#    cellnote=matrixAAin[,-takeout][,no],
 #   notecol="black",
 #   notecex=1.2,
#    adjCol=0.5,
	srtCol=45,
#    keysize=ks,key.title="",key.xlab="-log10(p-value)",key.ylab="",density.info="histogram",
    cexRow=1.7,cexCol=1.5,lwid=c(0.1,0.8))
dev.off()



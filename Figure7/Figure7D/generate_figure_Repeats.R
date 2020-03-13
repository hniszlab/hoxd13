#!/usr/bin/Rscript

# generate_figure_Repeats.R -  R script
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

#outmat=read.table("./data/outmatrix_v4.csv")
outmat=read.table("./outmatrix_v4.csv")
repenr<-table(outmat[,5])

## cluster sizes
csizes=table(outmat[,3])
reslist=list()
reslist.i=0
enrich.mat <-list()

## matrix counting the number of AA enrichments inside
matrixAAin=matrix(0,nrow=7,ncol=length(repenr))
colnames(matrixAAin)=names(repenr)



for(n in names(repenr)){
	reslist.i=reslist.i+1
	## the enrichment matrices
	enrich.mat[[reslist.i]]=list()
	resv=c()
	content_repenr <- table(as.numeric(outmat[which(outmat[,5]==n),3]))
	
	## foreach cluster we ask now 
	for(cin in 1:7){
		if(cin %in% names(content_repenr)){
			ci=which(names(content_repenr) == cin)
	
		cname=names(content_repenr)[ci]

		## index in csizes that matches our cluster
		clust.i <-which(names(csizes) == cname)

		## how many in current cluster
		camount=content_repenr[ci]
		## how many are in other clusters
		rest=sum(content_repenr)-camount

		## what is left in current cluster
		caleft=csizes[clust.i]-camount

		## what is left from total then is 

		total=nrow(outmat)-rest-camount-caleft


		mmat=matrix(c(camount,rest,caleft,total),ncol=2)
		colnames(mmat) =c("AAin","AAother")
		rownames(mmat) =c("ClutsterGnees","Othergenes")
		enrich.mat[[reslist.i]][[cin]]=mmat
		matrixAAin[cin,reslist.i] <-camount

	
			ft <- fisher.test(mmat,alternative="greater")
			resv=c(resv,ft[[1]])
		}else{
			resv=c(resv,1)
			enrich.mat[[reslist.i]][[cin]]=0
		}
	}
	reslist[[reslist.i]]=resv
	names(enrich.mat[[reslist.i]])=1:7
}
names(reslist)=names(repenr)
names(enrich.mat)=names(repenr)


resmat<-matrix(unlist(reslist),nrow=7)
colnames(resmat)=names(repenr)
rownames(resmat)=1:7
restmat_binary=resmat
restmat_binary[]=0
restmat_binary[resmat < 0.05] = 1


reord.disp <- c(1,6,3,7,5,4,2) ## this reoders the columns only !

## and the color gets now some pvalue color

colsin=foldc.colors9[5:9]
breaksin=c(0,1.30103,3,5,10,50)


mp <-function(of1=95,of2=125,svg=F){
	ks=1
	if(svg==T){
		ks=1.5
		par(oma=c(7,10,5,10))
	}
	heatmap.2(-log10(resmat[,reord.disp]),
			  dendrogram="none",Rowv=F,Colv=F,
			  tracecol=NA,
			  RowSideColors=colors7,srtCol=1,
			  #	xlab="Enrichment for members\ncontaining the indicated\nAA repeat",
			  offsetRow=-of1,offsetCol=-of2,
			  sepcolor="black",sepwidth=c(0.005,0.005),colsep=0:7,rowsep=0:7,
			  labRow=paste("IDR cluster",1:7),
			  #col=colorRampPalette(c("white","red"))(n = 8) ,
			  col=colsin,
			  #breaks=c(0,1.3,2,3,4,5,6,7,50),
			  breaks=breaksin,
			  labCol=paste(names(repenr),addbraces(repenr))[reord.disp],
			  cellnote=matrixAAin[,reord.disp],
			  notecol="black",
			  notecex=2,
			  adjCol=0.5,
			  keysize=ks,key.title="",key.xlab="-log10(p-value)",key.ylab="",density.info="histogram",
			  cexRow=2,cexCol=2)
}
svglite("Figure7D_repeats.svg",12,10.54687)
mp(svg=T,60,80);
mtext("Enrichment for members\ncontaining the indicated\nAA repeat",1,line=10,cex=2)
dev.off()

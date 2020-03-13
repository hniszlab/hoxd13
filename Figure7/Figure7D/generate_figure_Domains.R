#!/usr/bin/Rscript

# generate_figure_Domains.R -  R script
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
colsin=foldc.colors9[5:9]
breaksin=c(0,1.30103,3,5,10,50)
csizes=table(outmat[,3])

TFfam <- sort(table(outmat[,2]),decreasing=T)
TFenrich=list()
tf.i=0

matrixAAin=matrix(0,nrow=7,ncol=length(TFfam))
colnames(matrixAAin)=names(TFfam)


for(i in 1:7){
	## get the cluster ids
    ci=which(outmat[,3] == i)

    ## activator in cluster
	res.activ=c()
	for(tf.i in 1:length(TFfam)){ ## this list we iterate over
		actin=names(TFfam)[tf.i]         ## foreach element in list
	    aic<-which(outmat[ci,2] %in%actin)

	    ## activator not in cluster
		anic <- as.numeric(TFfam[tf.i])-length(aic)

		## index in csizes that is our cluster
		clust.i<-which(names(csizes) == i)
		## what is left in cluster that is not a Activ/Repress
		caleft=csizes[clust.i]-length(aic)

		## what is globally left
	    atotal=nrow(outmat)-length(aic)-anic-caleft
		matrixAAin[i,tf.i]=length(aic)

		mmat=matrix(c(length(aic),anic,caleft,atotal),ncol=2)
		ft <- fisher.test(mmat,alternative="greater")
	    res.activ=c(res.activ,ft[[1]])
	}
	TFenrich[[i]]=res.activ
}

resmat<-matrix(unlist(TFenrich),nrow=7,byrow=T)
colnames(resmat)=paste(names(TFfam),addbraces(TFfam))
rownames(resmat)=1:7
keep=c(3,4,2,1,6,7,8)

svglite("Figure7D_domains.svg",12,10)
par(oma=c(14,2,2,10))
ks=1
of1=0
of2=-10

heatmap.2(-log10(resmat[,keep]),
    dendrogram="none",Rowv=F,Colv=F,
    tracecol=NA,
    RowSideColors=colors7,
#   xlab="Enrichment for members\ncontaining the indicated\nAA repeat",
    offsetRow=-of1,offsetCol=-of2,
    sepcolor="black",sepwidth=c(0.005,0.005),colsep=0:7,rowsep=0:7,
    labRow=paste("IDR cluster",1:7),

    col=colsin,

    breaks=breaksin,
    cellnote=matrixAAin[,keep],
    notecol="black",
    notecex=2,
    adjCol=0.5,
    keysize=ks,key.title="",key.xlab="-log10(p-value)",key.ylab="",density.info="histogram",
    cexRow=2,cexCol=2,lwid=c(0.3,0.6))
dev.off()


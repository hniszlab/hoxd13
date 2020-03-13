#!/usr/bin/Rscript

# generate_figure_5D.R -  R script
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
DPATH="../../Fdata/Figure5D/"
SPATH="../../src/"
source(paste(SPATH,"colors.R",sep=""))
library(svglite)

library(RColorBrewer)
library(gplots)
a=as.matrix(read.table(paste(DPATH,"all_diffex_genes_adj05.tsv",sep=""),row.names=1))

a1=a
a1[a1>0]=1
a1[a1<0]=-1

minc=0
tt <- which(rowSums(abs(a)) > minc)
          #blue      greyish   "red"
mycols=c("#3954a5","#eff0f1","#ed2024")
mycols2=rev(mycols)

pdf("tmp123.pdf",13,13)
#hm<-heatmap.2(a1[rowSums(a1)!=0,],
hm<-heatmap.2(a1[tt,],
		  Colv=F,Rowv=T,dendrogram="none",trace="none",col=mycols2,labCol=1:11,main="diffex_genes_adjp_lt_0.05_expressed",
		  tracecol=NA,
		  density="density",
		  keysize=1,key.ylab=NA,
		  breaks=4,
		  srtCol=0,
		  sepcolor="",
         key.title=NA,
         key.xlab=NA,
         key.par=list(mgp=c(1.5, 0.5, 0),
                      mar=c(1, 2.5, 1, 0)),
         key.xtickfun=function() {
              cex <- par("cex")*par("cex.axis")
              side <- 1
              line <- 0
              col <- par("col.axis")
              font <- par("font.axis")
              mtext("Up SPDH", side=side, at=0, adj=0,
                    line=line, cex=cex, col=col, font=font)
              mtext("Down SPDH", side=side, at=1, adj=1,
                    line=line, cex=cex, col=col, font=font)
              return(list(labels=FALSE, tick=FALSE))
         })

		mtext("Cluster",1,7)
dev.off()
matt=a1[tt,]

ro <- order(rowSums(matt))
matt2=matt
matt2[matt2>0]=0
ro2 <- order(rowSums(matt2))
nmat=a1[tt,][hm$rowInd,]
omat=nmat

if(minc == 0){
	nmat[1,]=a1[tt,][hm$rowInd,][30,]
	nmat[2,]=a1[tt,][hm$rowInd,][31,]
	nmat[30,]=a1[tt,][hm$rowInd,][1,]
	nmat[31,]=a1[tt,][hm$rowInd,][2,]


	rownames(nmat)[1] = "Xist"
	rownames(nmat)[2] = "Rm45"
	rownames(nmat)[30] = rownames(omat)[1]
	rownames(nmat)[31] = rownames(omat)[2]
}

svglite("Figure5D.svg",13,20)

heatmap.2(nmat,
		  sepwidth=c(0,0),
		  sepcolor=NA,
          Colv=F,Rowv=F,dendrogram="none",trace="none",col=mycols2,labCol=1:11,main="diffex_genes_adjp_lt_0.05_expressed",
          tracecol=NA,
          density="density",
          keysize=1,key.ylab=NA,
          breaks=4,
          srtCol=0,
         key.title=NA,
         key.xlab=NA,
         key.par=list(mgp=c(1.5, 0.5, 0),
                      mar=c(1, 2.5, 1, 0)),
         key.xtickfun=function() {
              cex <- par("cex")*par("cex.axis")
              side <- 1
              line <- 0
              col <- par("col.axis")
              font <- par("font.axis")
              mtext("Up SPDH", side=side, at=0, adj=0,
                    line=line, cex=cex, col=col, font=font)
              mtext("Down SPDH", side=side, at=1, adj=1,
                    line=line, cex=cex, col=col, font=font)
              return(list(labels=FALSE, tick=FALSE))
         })

        mtext("Cluster",1,7)
dev.off()
write.table(rownames(nmat),paste("Figure5D_genes_in_heatmap_minclust_",(minc+1),".txt",sep=""),col.names=F,quote=F)

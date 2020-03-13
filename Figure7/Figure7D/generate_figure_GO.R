#!/usr/bin/Rscript

# generate_figure_GO.R -  R script
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
DPATH="../../data/"
SPATH="../../src/"

source(paste(SPATH,"colors.R",sep=""))

a=read.table("go_outtable.txt",row.names=1)
hin <-scan("go_outtable.txt",what="char",nlines=1)
colnames(a)=hin[-1]


#resmat=t(a[,8:14])     ## uncorrected pvalues
resmat=t(a[,1:7])      ## FDR corrected
matrixAAin=t(a[,22:28])
matrixAAin2=t(a[,15:21])

## number of total entries per category
tentries<-apply(a[,15:21],1,max)
colnames(resmat)=paste(colnames(resmat),addbraces(tentries))

takeout<-which(tentries < 20)                          ## filter on min set size
takeout<-c(takeout,which(tentries > 1000))             ## filter on max set size
takeout<-c(takeout,which(apply(resmat,2,min) >= 0.05)) ## filter on q/p -value
takeout<-c(takeout,which(apply(matrixAAin,2,max)< 15)) ## filter on abs number of enriched genes per set in cluster  

no <- order(apply(matrixAAin[,-takeout],2,max),decreasing=T)

svglite("Figure7D_GO.svg",15,10)
par(oma=c(30,20,1,5))
ks=1
of1=0
of2=0
heatmap.2(-log10(resmat[,-takeout][,no]),
    dendrogram="none",Rowv=F,Colv=F,
    tracecol=NA,
    offsetRow=-of1,offsetCol=-of2,
    labRow=paste("IDR cluster",1:7),
    col=foldc.colors9[5:9],
    breaks=c(0,1.30103,3,5,10,50),
   #cellnote=matrixAAin[,-takeout][,no],
   #notecol="black",
   #notecex=1.2,
    #adjCol=0.5,
    #srtCol=45,
    key.title="",key.xlab="-log10(p-value)",key.ylab="",density.info="histogram",
    cexRow=1.7,cexCol=1.5,lwid=c(0.1,0.8))

dev.off()


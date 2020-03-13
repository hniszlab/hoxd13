#!/usr/bin/Rscript

# generate_figure_Function.R -  R script
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

outmat=read.table(paste(DPATH,"Figure7/outmatrix_v4.csv",sep=""))

colsin=foldc.colors9[5:9]
breaksin=c(0,1.30103,3,5,10,50)

repressors<-scan(paste(DPATH,"Figure7/transcription_factors/repressors_uc.csv",sep=""),what="char")
activators<-scan(paste(DPATH,"Figure7/transcription_factors/activators_uc.csv",sep=""),what="char")

repin <- repressors[which(repressors %in% outmat[,1])]
actin <- activators[which(activators %in% outmat[,1])]

csizes=table(outmat[,3])



res.activ=c()
res.repres=c()
res.activ.l=list()
res.repres.l=list()
matrixAAin=matrix(0,nrow=7,ncol=2)

for(i in 1:7){
    ci=which(outmat[,3] == i)
    ## activator in cluster
    aic<-which(actin%in%outmat[ci,1])

    ## repressor in cluster
    ric<-which(repin%in%outmat[ci,1])

    ## repressor not in cluster
    rnic <- length(repin)-length(ric)
    ## activator not in cluster
    anic <- length(actin)-length(aic)

    ## index in csizes that is our cluster
    clust.i<-which(names(csizes) == i)
    ## what is left in cluster that is not a Activ/Repress
    caleft=csizes[clust.i]-length(aic)
    crleft=csizes[clust.i]-length(ric)

    ## what is globally left
    atotal=nrow(outmat)-length(aic)-anic-caleft
    rtotal=nrow(outmat)-length(ric)-rnic-crleft

    mmat=matrix(c(length(aic),anic,caleft,atotal),ncol=2)
    mmatr=matrix(c(length(ric),rnic,crleft,rtotal),ncol=2)

    ft <- fisher.test(mmat,alternative="greater")
    ftl <- fisher.test(mmat,alternative="less")
    #ft <- fisher.test(mmat)
    res.activ=c(res.activ,ft[[1]])
    res.activ.l<-c(res.activ.l,ftl[[1]]) ## add depleted
    ft <- fisher.test(mmatr,alternative="greater")
    ftl <- fisher.test(mmatr,alternative="less")
    #ft <- fisher.test(mmatr)
    res.repres=c(res.repres,ft[[1]])
    res.repres.l<-c(res.repres.l,ftl[[1]]) ## add depleted
    matrixAAin[i,1:2]=c(length(aic),length(ric))
}


resmat=cbind(res.activ,res.repres)
colnames(resmat)=c(paste("Activators",addbraces(length(actin))),paste("Repressors",addbraces(length(repin))))

svg("Figure7D_Function.svg",7,10)
par(oma=c(12,5,4,12))

ks=1
of1=0
of2=-10
heatmap.2(-log10(resmat),
    dendrogram="none",Rowv=F,Colv=F,
    tracecol=NA,
    RowSideColors=colors7,
#   xlab="Enrichment for members\ncontaining the indicated\nAA repeat",
    offsetRow=-of1,offsetCol=-of2,
    sepcolor="black",sepwidth=c(0.005,0.005),colsep=0:7,rowsep=0:7,
    labRow=paste("IDR cluster",1:7),
    #col=colorRampPalette(c("white","red"))(n = 8) ,
    col=colsin,
    #breaks=c(0,1.3,2,3,4,5,6,7,50),
    breaks=breaksin,
    cellnote=matrixAAin,
    notecol="black",
    notecex=2,
    adjCol=0.5,
    keysize=ks,key.title="",key.xlab="-log10(p-value)",key.ylab="",density.info="histogram",
    cexRow=2,cexCol=2,lwid=c(0.4,0.6))
dev.off()


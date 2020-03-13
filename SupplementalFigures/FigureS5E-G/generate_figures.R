#!/usr/bin/Rscript

# generate_figure_S5E-G.R -  R script
# Copyright (C) 2018 Sebastian Mackowiak
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
DPATH="../../Fdata/FigureS5E-G/"
SPATH="../../src/"


library(svglite)

## read in file with hoxd13, gene, TAD
all4 <-read.table(paste(DPATH,"TADS_with_hoxd13_cluster4",sep=""))

#now load all TADs with ctrl genes
allo <-read.table(paste(DPATH,"tads_with_genes_ctrl",sep=""))

a2=allo[which(allo[,2] >0),]

## FigureS5E
a=c(80/82,2755/3563)
svg("FigureS5E.svg",6,12)
barplot(a,las=1,yaxt="none",col=c("brown","darkblue"),
        names.arg=c("Cluster4 genes","Non-cluster4\ngenes"),
        ylab="ratio of TADS with gene and hoxd13 peak",
        ylim=c(0,1))
axis(2,seq(0,1,0.1),paste(seq(0,100,10),"%",sep=""),las=1,ylab="ratio of TADS with gene and hoxd13 peak")
dev.off()


# Figure S5F
resn=c();
for(i in 1:1000){
        od=sample(1:dim(allo)[1]);
        resn=c(resn,length(which(allo[od[1:82],2] > 0)))
}

svglite("FigureS5F.svg",5,10)
plot(density(resn/82),xlim=c(0.5,1),xlab="ratio [Tads with hoxd13/all Tads]",main="sampling 82 tads a 1000 times from all tads\nthe empirical p-value is 0 since none of the 1000 ratios from shuffling\n was as good as the real one from cluster 4 (blue line)",yaxt="n",ylim=c(0,10),col="orange",ylab="Frequency")
abline(v=80/82,col="#b7e7ff",lwd=2)
axis(2,seq(0,10,2),seq(0,500,100),las=1)
dev.off()


## FigureS5G
dista <- read.table(paste(DPATH,"closest_hox_allo.csv",sep=""))
dist4=as.matrix(read.table(paste(DPATH,"closest_hox_with_tad_etc.csv",sep="")))
resX=c()
for(i in 1:1000){
	od=sample(1:dim(dista)[1])
	resX=c(resX,dista[od[1:84],2])
}


svglite("FigureS5G.svg",20,7)
ylims1=c(0,14)
ylims2=c(0,2e4)
yl="Number of occurrence"

par(mfrow=c(1,2))
xran <- c(0,24)
hist(log2(as.numeric(dist4[,4])),
	 freq=T,breaks=10,
	 col="cornflowerblue",las=1,
	 xlab="Hoxd13 peak distance from TSS [nt]",
	 main="Distance from the TSS to the closest Hoxd13 peak\nCluster 4 dysregulated genes (n=84)",
	 xaxt="n",
	 ylim=ylims1,
	 xlim=xran,
	 ylab=yl)

ran=seq(0,24,2)
ggo<-c(2^ran)
ggg<-round(ggo/1024)
xlabx=c(ggo[which(ggg==0)], paste(ggg[which(ggg < 1000 & ggg > 0)],"k",sep=""),paste(round(ggg[which(ggg > 1000)])/1024,"M",sep=""))
axis(1,ran,xlabx)
abline(v=log2(median(as.numeric(dist4[,4]))),col="brown",lwd=2)


hist(log2(resX),
	 freq=T,breaks=20,
	 col="brown",las=1,
	 xlab="Hoxd13 peak distance from TSS [nt]",
	 main="Distance from the TSS to the closest Hoxd13 peak\n84 randomly selected genes (iterated 1000 times)",
	 xaxt="n",
	 ylim=ylims2,
	 xlim=xran,
	 ylab=yl)

axis(1,ran,xlabx)
d2<-density(log2(resX))
to2<-which(d2$x < 1)
abline(v=log2(median(resX)),col="#b7e7ff",lwd=2)
dev.off()

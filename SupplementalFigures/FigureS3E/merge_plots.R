#!/usr/bin/Rscript

# merge_plot.R -  R script
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


library(raster)
library(RColorBrewer)
library(scales)
library(data.table)

## Helper functions
## join strings
strjoin <- function(x,sep="_"){
	paste0(x,collapse=sep)
}

## metaplot
metaPlot <- function(mat,...){
	par(cex=1.3,mar=c(mar[1:3],0.5))
	yv=colMeans(mat)
	plot(1:ncol(mat),yv,"l",ylab="average rpm/bp",lwd=2,col="cornflowerblue",las=1,xlab="",xaxt="n",xaxs="i",yaxs="i",bty="n",...)
	axis(1,c(1,50,100,150),c(lflank,"S","E",rflank))
}

### HEATMAP ###
heatPlot <- function(mat,axy=T,mplot=T,allF=F,...){
	if(!allF){
		par(cex=1.3,mar=c(mar[1:2],0.5,0.5))
	}

	ylab=""
	if(axy){
		ylab="Peaks"
	}
	#image(t(fmat2),las=1,xaxt="n",yaxt="n",ylab=ylab,bty="n",...)
	if(mplot){
		raster::image(t(mat),las=1,xaxt="n",yaxt="n",ylab=ylab,bty="n",...)
	}else{
		nplot(xlim=c(0,1),ylim=c(0,1))
	}
	if(!allF){
		axis(1,c(0,0.333,0.666,1),c(lflank,"S","E",rflank),lwd=0,lwd.ticks=1)
	}
	if(axy){
		axis(2,seq(0,1,length.out=11),ceiling(seq(nrow(mat),1,length.out=11)),las=1,lwd=0,lwd.ticks=1)
	}
}

### color bar
colBar <- function(mat,colUse=hcl,rm=3,mplot=T,...){
	ymi=min(mat)
	yma=max(mat)
	grid <- raster(ncols=1, nrows = 1000, xmn=0, xmx=1, ymn=ymi,ymx=yma,vals=seq(ymi,yma,length.out=1000))
	
	## the color bar on the side
	par(mar=c(mar[1],0,0.5,rm))
	nplot(0,xlim=c(0,1),ylim=c(ymi,yma),xaxs="i",yaxs="i")
	yr=seq(ymi,yma,length.out=10)
	axis(4,yr,round(yr,2),las=1,lwd=0,lwd.ticks=1)
	plot(grid, col=rev(colUse), legend=FALSE, axes = 0, box=F,nc=1,nr=1000,add=T)
	mtext("rpm/bp",side=4,line=1,cex=1.5,at=par("usr")[3]+0.4*diff(par("usr")[3:4]),outer=T)
}
##################

args = commandArgs(trailingOnly=TRUE)

name1=strsplit(args[1],"_")[[1]]
name2=strsplit(args[2],"_")[[1]]
pref=args[3]

ln1=length(name1)
ln2=length(name2)

flank=as.numeric(args[4])
flank=flank/1000

if(length(args) < 5){
	args[5]="BuRd"
}



## they are standardized 
sufMeta='intersect_tag_fused_sorted_meta_plot.csv'
sufHeat='intersect_tag_fused_sorted_heatmap.csv'

p1=paste0(args[1],"/",args[3])
p2=paste0(args[2],"/",args[3])

ids=scan(paste0(args[1],"/","regions_ids"),what="char",skip=1)
fused=grep("fused",ids)
S1=grep(name1[(ln2-1)],ids)
S2=(1:length(ids))[-c(fused,S1)]



#f1meta=as.matrix(read.table(paste(p1,sufMeta,sep=".")))
#f2meta=as.matrix(read.table(paste(p2,sufMeta,sep=".")))
f1meta=fread(paste(p1,sufMeta,sep="."))
f2meta=fread(paste(p2,sufMeta,sep="."))

quants=as.numeric(quantile(unlist(f1meta),c(0.01,0.98)))

p01<-quants[1]
p99<-quants[2]

ss=1:nrow(f1meta)

## just using the fused peaks as default, when using all regions we still get a signal of the
## other sample which has no detected peak ( does not mean nothing was detected, simply no peak was called)
## thus we see signal in both samples
all=ss
for(tsel in c("all")){ #,"fused","S1","S2")){
#for(tsel in c("fused")){

ss=get_var(tsel)
outadd=paste0("_",tsel)


## use a subset here now which is in ss
fmat1=as.matrix(f1meta[ss,])
no1=order(rowMeans(fmat1[,61:90]))
fmat1[fmat1>p99]=p99
fmat1[fmat1<p01]=p01

quants2=as.numeric(quantile(unlist(f2meta),c(0.01,0.98)))
p01<-quants[1]
p99<-quants[2]

fmat2=as.matrix(f2meta[ss,])
no2=order(rowMeans(fmat2[,61:90]))
fmat2[fmat2>p99]=p99
fmat2[fmat2<p01]=p01

f1heat=fmat1[no1,]
## use no2 to get original ordering 
f2heat=fmat2[no1,] ## use same order as we have for the first matrix, so we can directly compare peaks in the heatmap

RdBu = brewer.pal(11,'RdBu');
BuRd = rev(RdBu)

bluered<-colorRampPalette(BuRd)
redblue<-colorRampPalette(RdBu)
viridis=viridis_pal()

un="kb"
if(flank > 1000){
    flank=flank/1000
    un="mb"
}

lflank=paste("-",flank,un,sep="")
rflank=paste("+",flank,un,sep="")

breaks=seq(min(f1heat),max(f1heat),length.out=1001)
#breaks2=seq(min(f1heat),2,length.out=1001)

## adjust image height
pixH=1400
pixH2=3
pixL2=1
if(nrow(f1heat) > 40000){
#	pixH=6000
#	pixH2=8
#	pixL2=1
}

colUse=bluered(1000)
if(args[5] == "viridis"){
    colUse=viridis(1000)
}

if(args[5] == "RdBu"){
    colUse=redblue(1000)
}

mar=c(5,6,4,2)+0.1
usecol=c(25:125)


png(paste(args[3],outadd,"S1_heatmap_",args[5],".png",sep=""),5,10,units="in",res=600)
par(mar=c(0,0,0,0))
heatPlot(f1heat[,usecol],col=colUse,axy=F,breaks=breaks,mplot=T,allF=T)

dev.off()

png(paste(args[3],outadd,"S2_heatmap_",args[5],".png",sep=""),5,10,units="in",res=600)
par(mar=c(0,0,0,0))
heatPlot(f2heat[,usecol],col=colUse,axy=F,breaks=breaks,mplot=T,allF=T)
dev.off()

#png(paste(args[3],outadd,".png",sep=""),800,pixH)
for(mplot in c(T,F)){
	if(mplot){
		png(paste(args[3],outadd,"_",args[5],".png",sep=""),10,14,units="in",res=300)
	}else{
		pdf(paste(args[3],outadd,"_",args[5],".pdf",sep=""),10,14,useDingbats=F)
	}
par(oma=c(0,0,3,2),mgp=c(4,1,0),mar=mar)
lmat=matrix(c(1,2,3,4,0,5),ncol=3,nrow=2)
layout(lmat,widths=c(4,4,1),heights=c(pixL2,pixH2))

ym1=colMeans(f1meta[ss,])
ym2=colMeans(f2meta[ss,])
ymm=max(ym1,ym2)
ymD=1.2
if(ymm > ymD){
	ymD=ymm
}

metaPlot(f1meta[ss,],ylim=c(0,ymD),main=name1[ln1])
lines(1:ncol(f2meta),colMeans(f2meta[ss,]),col="grey",lwd=2)
heatPlot(f1heat[,usecol],col=colUse,mplot=mplot)

metaPlot(f2meta[ss,],ylim=c(0,ymD),main=name2[ln2])
lines(1:ncol(f1meta),colMeans(f1meta[ss,]),col="grey",lwd=2)
heatPlot(f2heat[,usecol],col=colUse,axy=F,breaks=breaks,mplot)

colBar(f1heat,colUse=colUse,rm=3,mplot=mplot)

## add common title
mtext(paste(args[3]," ",tsel," peaks\nFlanks = ",flank,un,sep=""),3,outer=T,cex=2)
dev.off()
}
}

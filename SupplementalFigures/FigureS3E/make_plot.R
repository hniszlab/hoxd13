#!/usr/bin/Rscript

# make_plot.R -  R script
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

library(argparse)
library(data.table)

## parse arguments first
parser <- ArgumentParser()
parser$add_argument("-r", "--rpmbp", default="NA",type="character",
    help="Filename of one of the output.tsv files from the rpmbp pipeline")

parser$add_argument("-b", "--bins", default=10,type="integer",
    help="Number of bins to remove from the front and end of each matrix. These are the overlapping ones")

parser$add_argument("-n", "--name", default="outmatrix",type="character",
    help="Name of the output files")

parser$add_argument("-t", "--title", default="Meta_and_Heatmap",type="character",
    help="plot_title")

parser$add_argument("-c", "--color", default="RdBu",type="character",
    help="Color scheme to be used. Currently implemented are RdBu,BuRd and viridis. Adding more later on.")

parser$add_argument("-f", "--flank", default=3000,type="integer",
    help="Size of the flanking regions in nt")

parser$add_argument("-z", "--nfac", default=1,type="double",
    help="Normalize values by factor")

args <- parser$parse_args()


## give one of the three files
if(args$rpmbp == "NA"){
	print("No file given for reading in")
	quit()
}


## get filename without extension
file=sub('\\..sv$', '', args$rpmbp)

library(raster)
library(RColorBrewer)
library(scales)

## determine file structure
check=c("_left_","_center_","_right_")
ss=c()
for(i in check){
	if(length(grep(i,file))>0){
		ss=strsplit(file,i)
	}
}

rem=10   ## this should not be changed unless you know what you are doing. If you adapt this you need to adapt the flank region creation for the center bins as well
if(args$bins){
	rem=as.numeric(args$bins)
}

### read in tables
l=list()
j=0
for(i in check){
	j=j+1
	l[[j]]=read.table(paste(ss[[1]][1],i,ss[[1]][2],".tsv",sep=""),skip=1)
}

regionsIds=l[[2]][,1]

names(l)=check
## this is 52 if we used 50 bins initially, first two columns are feature descriptions
bins=ncol(l[[2]])

## in meta_plots.sh
br=3+rem
be=bins-rem
lmat=l[[1]][,br:be]
cmat=l[[2]][,br:be]
rmat=l[[3]][,br:be]

## now exclude peaks that do not have full flanks covered
ex=c()
ex=c(ex,unique(which(lmat == -1, arr.ind=T)[,1]))
ex=c(ex,unique(which(rmat == -1, arr.ind=T)[,1]))
ex=c(ex,unique(which(cmat == -1, arr.ind=T)[,1]))
ex=unique(ex)

if(length(ex) > 0){
lmat=lmat[-ex,]
rmat=rmat[-ex,]
cmat=cmat[-ex,]
regionsIds=regionsIds[-ex]
}
write.table(regionsIds,file="regions_ids",quote=F,row.names=F,col.names=F)

fmat=as.matrix(cbind(lmat,cmat,rmat))

mm=ss[[1]][2]
plotTitle=args$title
flank=args$flank/1000

un="kb"
if(flank > 1000){
	flank=flank/1000
	un="mb"
}

lflank=paste("-",flank,un,sep="")
rflank=paste("+",flank,un,sep="")


yl=c(0,max(colMeans(fmat)+0.1))
pdf(paste(mm,"_meta.pdf",sep=""),10,10,useDingbats=F)
plot(1:ncol(fmat),colMeans(fmat),"l",ylab="average rpm/bp",lwd=2,col="cornflowerblue",las=1,xlab="bins",main=mm,ylim=yl,xaxt="n")
axis(1,c(1,50,100,150),c(lflank,"S","E",rflank))
dev.off()

write.table(fmat,file=paste(mm,"_meta_plot.csv",sep=""),sep="\t",quote=F,row.names=F,col.names=F)

## adjust this value to get more or less background color
p99<-as.numeric(quantile(unlist(fmat),0.98))
p01<-as.numeric(quantile(unlist(fmat),0.01))

RdBu = brewer.pal(11,'RdBu');
BuRd = rev(RdBu)

bluered<-colorRampPalette(BuRd)
redblue<-colorRampPalette(RdBu)
viridis=viridis_pal()

colUse=bluered(1000)
if(args$color == "viridis"){
	colUse=viridis(1000)
}

if(args$color == "BuRd"){
	colUse=redblue(1000)
}



## make the heatmap
#no=order(rowMeans(fmat[,51:100]))
no=order(rowMeans(fmat[,61:90]))
fmat2=as.matrix(fmat)
fmat2[fmat2>p99]=p99
fmat2[fmat2<p01]=p01
png(paste(mm,"_heatmap.png",sep=""),200,1024)
#image(t(fmat2[no,]),col=bluered(1000),las=1,xaxt="n")
image(t(fmat2[no,]),col=colUse,las=1,xaxt="n")
axis(1,c(0,0.333,0.666,1),c(lflank,"S","E",rflank))
dev.off()

write.table(fmat2[no,],file=paste(mm,"_heatmap.csv",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
write.table(no,file=paste(mm,"_columnOrder.csv",sep=""),sep="\t",quote=F,row.names=F,col.names=F)

## orig mar
mar=c(5,5,4,2)+0.1

#png(paste(mm,"_meta_heatmap.png",sep=""),500,1400)
png(paste(args$title,"_meta_heatmap.png",sep=""),8,14,units="in",res=300)

## define plot matrix, plotting is column wise, and the upper right corner is not used in here
lmat=matrix(c(1,2,0,3),2,2)
layout(lmat,widths=c(4,1),heights=c(1.3,3))

### META ### 
par(cex=1.3,mar=c(mar[1:3],0.5))
plot(1:ncol(fmat),colMeans(fmat),"l",ylab="average rpm/bp",lwd=5,col="cornflowerblue",las=1,xlab="",main=plotTitle,ylim=yl,xaxt="n",xaxs="i",yaxs="i",bty="n")
axis(1,c(1,50,100,150),c(lflank,"S","E",rflank))

### HEATMAP ###
par(cex=1.3,mar=c(mar[1:2],0.5,0.5),mgp=c(4,1,0))
image(t(fmat2[no,]),col=colUse,las=1,xaxt="n",yaxt="n",ylab="Peaks",bty="n")
axis(1,c(0,0.333,0.666,1),c(lflank,"S","E",rflank),lwd=0,lwd.ticks=1)
axis(2,seq(0,1,length.out=11),ceiling(seq(nrow(fmat2),1,length.out=11)),las=1,lwd=0,lwd.ticks=1)

### COLOR_BAR ###
ymi=min(fmat2)
yma=max(fmat2)
grid <- raster(ncols=1, nrows = 1000, xmn=0, xmx=1, ymn=ymi,ymx=yma,vals=seq(ymi,yma,length.out=1000))
## the color bar on the side
par(mar=c(mar[1],0,0.5,5))
nplot(0,xlim=c(0,1),ylim=c(ymi,yma),xaxs="i",yaxs="i")
yr=seq(ymi,yma,length.out=10)
axis(4,yr,round(yr,2),las=1,lwd=0,lwd.ticks=1)
plot(grid, col=rev(colUse), legend=FALSE, axes = 0, box=F,nc=1,nr=1000,add=T)
mtext("rpm/bp",side=4,line=3,cex=1.5,at=par("usr")[3]+0.5*diff(par("usr")[3:4]))
dev.off()

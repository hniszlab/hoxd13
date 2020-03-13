#!/usr/local/bin/Rscript

# make_raster_plot.R -  R script
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

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
	args=c("rpmbp_wt1_vs_wt2.intersect_tag_fused_sorted.bed_wt1","rpmbp_wt1_vs_wt2.intersect_tag_fused_sorted.bed_wt2","wt1","wt2")
}

sample1=read.table(args[1],header=T)
sample2=read.table(args[2],header=T)

x=sample1[,3]
y=sample2[,3]

ct <- cor.test(log2(x),log2(y),method="spearman")
		   
p.mar<-par("mar")
dev.off()
## dont print exponentials
options("scipen"=100, "digits"=4)

in1="\n(log10[rpm per bp])"
in2="\n(RPM per bp)"

inf=in2
if(0){
pdf(paste(args[3],"_vs_",args[4],".pdf",sep="x"))
par(mar=c(5.1,5.1,2.1,2.1))
plot(log10(x),log10(y),pch=16,col=adjustcolor( "blue", alpha.f = 0.1),
	 xlab=paste(args[3],inf),
	 ylab=paste(args[4],inf),
	 las=1,
	 xlim=c(-3,3),
	 ylim=c(-3,3),
	 axes=F
)

axv <- seq(-3,3,1)
axvec<-signif(10^axv,2)

axis(1,axv,axvec)
axis(2,axv,axvec,las=1)
abline(0,1,col="grey",lty=2)
text(2,-2.8,paste("Rho =",as.numeric(round(ct[[4]],2))),col="red",cex=2)
dev.off()

quit()
}

library(png)
#pdf("test.pdf", width = 7, height = 7)
pdf(paste(args[3],"_vs_",args[4],"raster.pdf",sep=""))
par(mar=c(5.1,5.1,5.1,5.1))
xl<-c(-3,3)
yl<-c(-3,3)

if(length(args) > 5){
	xl=c(as.numeric(args[5]),as.numeric(args[6]))
	yl=xl
}

plot(log10(x),log10(y),pch=16,col=adjustcolor( "blue", alpha.f = 0.1),
     xlab=paste(args[3],inf),
     ylab=paste(args[4],inf),
     las=1,
     xlim=xl,
     ylim=yl,
     axes=F,
	 type="n"
)
axv <- seq(xl[1],xl[2],1)
axvec<-signif(10^axv,2)

axis(1,axv,axvec)
axis(2,axv,axvec,las=1)
#abline(0,1,col="grey",lty=2)
if(xl[1] == -3){
	text(2,-2.8,paste("Rho =",as.numeric(round(ct[[4]],2))),col="red",cex=2)
}else{
	text(0.5,-0.8,paste("Rho =",as.numeric(round(ct[[4]],2))),col="red",cex=2)
}
# Extract plot area in both user and physical coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

tmp <- tempfile()

png(tmp, width = width, height = height, units = "in", res = 600, bg = "transparent")
#png(tmp, width = width, height = height, units = "in", res = 300)
par(mar=c(0,0,0,0))
plot(log10(x),log10(y),pch=16,col=adjustcolor( "blue", alpha.f = 0.1),axes=F,xlim=xl,ylim=yl)
abline(0,1,col="grey",lty=2)
dev.off()

# Windows users may have trouble with transparent plot backgrounds; if this is the case,
# set bg = "white" above and move the legend plot command below the raster plot command.
panel <- readPNG(tmp)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

dev.off()

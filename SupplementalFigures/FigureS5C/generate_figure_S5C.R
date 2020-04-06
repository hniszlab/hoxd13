#!/usr/bin/Rscript

# generate_figure_S5C -  R script
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


nums=read.table("overlaps_cluster9.nums")
library(fmsb)
## if the installation of the package by install.packages fails then 
## download the package from CRAN https://cran.r-project.org/web/packages/fmsb/index.html
## and install it manually by 
## R CMD INSTALL fmsb_0.7.0.tar.gz 
library(svglite)

svglite("Figure_S5C.svg",14,7)
#png("Figure_S5C.png",600,600)
par(cex=1.6)

#myd2<-matrix(c(rep(0.45,11),rep(0,11),nums[,6]),nrow=3,byrow=T)
myd2<-matrix(c(rep(0.9,11),rep(0,11),nums[,6]),nrow=3,byrow=T)
colnames(myd2)=paste("cluster",1:11,sep="")
data2=data.frame(myd2)
data=data2[,c(1,11:2)]


radarchart( data  , axistype=1 ,
    #custom polygon
	pfcol=rgb(127, 176, 255,200,maxColorValue=255) ,
    pcol=rgb(127, 176, 255,maxColorValue=255) , plwd=1,pty="",

    #custom the grid
    cglcol="grey", cglty=1, axislabcol="darkblue", caxislabels=paste(seq(0,90,10),"%"), cglwd=0.8,
    seg=9,
    #custom labels
    vlcex=1,
	title("overlap of genes from clusters 1 to 11\nwith genes from cluster 9") ## it's cluster 9 counting from 1, original 8 counting from 0
)
dev.off()

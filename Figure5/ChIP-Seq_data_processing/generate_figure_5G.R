#!/usr/bin/Rscript

# generate_figure_5G.R -  R script
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
DPATH="../../Fdata/Figure5G/"
SPATH="../../src/"
source(paste(SPATH,"colors.R",sep=""))
library(svglite)
#a=c(80/82,2755/3563)

## read in file with hoxd13, gene, TAD
all4 <-read.table(paste(DPATH,"TADS_with_hoxd13_cluster4",sep=""))

#now load all TADs with ctrl genes
allo <-read.table(paste(DPATH,"tads_with_genes_ctrl",sep=""))
a2=allo[which(allo[,2] >0),]

## shuffle a 1000 times 82 TADs which have hoxd13 and TSS and take mean value of this shuffle
resA=c()
for(i in 1:1000){
	od=sample(1:dim(a2)[1])
	resA=c(resA,mean(a2[od[1:82],2]))
}

## Figure5G
svglite("Figure5G.svg",4,12)
par(mar=c(8,5,2,1),bty="n",cex=1.2)
nl=list(all4[,2],resA)
names(nl)=c("C4","Ctrl")
boxplot(nl,las=1,xlab="",
		main="Average Hoxd13 peaks\n per TAD",
		ylab="Hoxd13 peaks per TAD",
		col=c("cornflowerblue","brown"),
		bty="n",
		ylim=c(0,50))
dev.off()


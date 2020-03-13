#!/usr/bin/Rscript

# generate_figure_5H.R -  R script
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
DPATH="../../Fdata/Figure5H/"
SPATH="../../src/"
source(paste(SPATH,"colors.R",sep=""))
library(svglite)

HL<-read.table(paste(DPATH,"HL_testout.R",sep=""))
MB<-read.table(paste(DPATH,"MB_testout.R",sep=""))

svg("Figure5H.svg",12,6)
par(mfrow=c(1,2))
mycolors=c("#a65628","#c49a6c","#add8e6","#00008b","#377eb8","#ff7f00","#fbb040","#f9ed32","#4daf4a","#f881bf","#e41a1c")
add="- total signal normalized by mean"
mm="\nhoxd13 summit at center of bin 11"
barplot(HL[,2]/mean(HL[,2]),main=paste("HL",add,mm),ylim=c(0,1.3),col=mycolors[2],las=1,names.arg=HL[,1],xlab="bin (5kb each)")
barplot(MB[,2]/mean(MB[,2]),main=paste("MB",add,mm),ylim=c(0,1.3),col=mycolors[3],las=1,names.arg=HL[,1],xlab="bin (5kb each)")
dev.off()


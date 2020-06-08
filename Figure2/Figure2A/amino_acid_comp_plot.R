#!/usr/local/bin/Rscript

# amino_acid_comp_plot.R -  R script
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


args <- commandArgs(TRUE)
#args=c("HOXD13.Rscript","A","HOXD13_emboss.R")
source(args[1])

ax<-seq(50,ceiling(slen/10)*10,50)
if(max(ax) < slen){
  ax=c(ax,max(ax)+50)
}

## make list with aa to highlight
aaAi=0
if(length(args) > 1){
  aaAi=which(aa == args[2])
}

plot(1:max(ax),col="white",ylim=c(0,20),xlab="Position",ylab="Amino acids",main="",bty="o",axes=F,xlim=c(1,max(ax)))
## y axis
axis(2,(1:20)-0.5,aa,las=2,tick=F)
axis(1,c(1,ax),c(),las=1)


ms=0.9
ms1=1-ms

for(i in 1:20){
  colA="black"
  if(i %in% aaAi){ colA="brown"}
  tp=sl[[1]][[1]][[i]]
  if(tp[1] > 0){
    arrows(tp,i-ms,tp,i-ms1,length=0,lwd=2,col=colA)
  }else{
    #arrows(tp,i-1,tp,i,length=0,lwd=2,col="green")
  }
}

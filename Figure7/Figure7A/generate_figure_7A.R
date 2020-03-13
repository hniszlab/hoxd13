#!/usr/local/bin/Rscript

# generate_figure_7A.R -  R script
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
DPATH="../../Fdata/"
SPATH="../../src/"

library(dendextend)
library(svglite)
source(paste(SPATH,"colors.R",sep=""))
source(paste(SPATH,"functions.R",sep=""))
source("labbook_analysis_functions.R")
source("circplot_fun_re.R_new")

## adapt so it an be easily used by all 
pos_set=read.table(paste(DPATH,"Figure7/lambert/ENSG_TRUE_tf_lambert_zf_split_AA_repeat_lengths_8plus.csv",sep=""))
pos_set10=read.table(paste(DPATH,"Figure7/lambert/ENSG_TRUE_tf_lambert_zf_split_AA_repeat_lengths_10plus.csv",sep=""))

DBD=paste(DPATH,"Figure7/lambert/2018_AddAlignments/",sep="")

## all sets now length filtered
tf.list  <- readAAm(sel=1:32,insets=tfsets)
dbd.list <- readAAm(sel=1:4,insets=dbdsets,P=DBD,addf='')

pmat.tf2_10 <- process.mat(tf.list,"sAA",thresAA=0.10)
pmat.dbd_10 <- process.mat(dbd.list,"sAA",thresAA=0.10) ## we have many less TF left with single AA enrichment for DBDs
save(pmat.tf2_10,file=paste(DPATH,"Figure7/pmat.tf2_10.RData",sep=""))
save(pmat.dbd_10,file=paste(DPATH,"Figure7/pmat.dbd_10.RData",sep=""))

tlax<-circplot(pmat=pmat.tf2_10,tfl=tf.list,pos_set=pos_set10,nfam=8,nclust=7,dims.use=10)

## figure7A
#png("Figure7A.png",1200,1200)
svglite("Figure7A.svg",12,12)
new_gene_order<-replot_circle(tlax,plotdend=T,cexl=1.4,colsFam=j.colors10,colsClust=colors7[c(1,4,3,5,6,2,7)])
dev.off()



## make outmatrix table ... 
ids = paste(tf.list$fullmatrix[,2],tf.list$fullmatrix[,3],sep="_")
allind=pmat.tf2_10$indices
pset=tlax[[6]][,7]
repeats<-tf.list$fullmatrix[pmat.tf2_10$indices,c(1,33:34)]
rep10 = which(as.numeric(repeats[,3]) >= 10)

outmat=matrix(nrow=1446,ncol=7)
outmat[,1]=rownames(pmat.tf2_10[[1]])
outmat[,2]=tf.list$sets[tlax[[5]]]
outmat[,3]=new_gene_order
outmat[pset,4]= as.vector(tlax[[6]][,6]) 
outmat[rep10,5]=repeats[rep10,2]
outmat[,6]=ids[allind]
outmat[,7]=tf.list$fullmatrix[allind,32]

outmat.ud <- cbind(outmat,tlax[[4]])

colnames(outmat.ud)=c("TF","DBD","cluster","AArepeat_TF","AA_repeat_IDR","ID","sequence","clusters.orig")
write.table(outmat.ud,"outmatrix_v4.csv",col.names=T,quote=F)







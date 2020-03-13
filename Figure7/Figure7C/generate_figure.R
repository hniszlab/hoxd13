#!/usr/bin/Rscript

# generate_figure.R -  R script
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

DPATH="../../Fdata/"
SPATH="../../src/"



ttype=c("IDR","DBD")
chromlist=list()
j=0
for(tt in ttype){
	j=j+1
	chromlist[[j]]=list()
	if(tt == "IDR"){
		print("reading IDR now")
		#b=read.table("/project/hnisz_lab_productive/IDRs/all_hg19_proteins/outmatrix_v4.csv_Rtable.RZF")
		b=read.table(paste(DPATH,"outmatrix_v4.csv_Rtable.R",sep=""))
		refseq=4
	}else{
		print("reading DBD now")
		## this is for the DBDs only now
		b=read.table(paste(DPATH,"all_DBD_ids.txt_Rexons.R",sep=""))
		refseq=7
	}

	#for(chr in 1:22)
	for(chr in c(1:22,"X","Y","M")){ ## added on 20190916
		
		sel=which(b[,1] == paste("chr",chr,sep=""))
		a=fread(paste("chr",chr,".phyloP100way.wigFix_ra",sep=""))
		names(a)="v"

		outv=c()
		genelist=list()
		nnames=c()
		if(length(sel)>0){
		for(i in 1:length(sel)){
			if(i==1){
				numl=1
				gene=as.character(b[sel[i],refseq])
				nnames[numl]=gene
			}
			outv=c()
			beg=b[sel[i],2]
			end=b[sel[i],3]
			outv=a[beg:end,as.numeric(v)]
			if(gene == as.character(b[sel[i],refseq])){
				if(i==1){
					genelist[[numl]]=c(outv)
				}else{
					genelist[[numl]]=c(genelist[[numl]],outv)
				}
			}else{
				gene=as.character(b[sel[i],refseq])
				numl=numl+1
				nnames[numl]=gene
				genelist[[numl]]=c(outv)
			}
		}
		}
		names(genelist)=nnames

		## get mean values
		mymeans <- lapply(genelist,mean)
		chromlist[[j]][[chr]]=unlist(mymeans)
	}
}

names(chromlist)=tt

chromlistIDR=chromlist[[1]]
chromlistDBD=chromlist[[2]]

ol=list(unlist(chromlistDBD),unlist(chromlistIDR))
names(ol)=c("DBD","IDR")

svglite("Figure7C.svg",3,6)
boxplot(ol,ylab="phyloP score",col=c("grey","steelblue4"),bty="n",las=1,outline=F)
dev.off()

write.table(ol[[1]],"phyloP_DBD.tsv",quote=F,col.names=F)
write.table(ol[[2]],"phyloP_IDR.tsv",quote=F,col.names=F)



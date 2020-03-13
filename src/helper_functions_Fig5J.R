#!/usr/local/bin/Rscript


# helper_functions_Fig5J.R -  R script
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


## take a list object and return a data frame object
listToDataFrame <- function(l){
	lenl=length(l)
	df=l[[1]]
	for(i in 2:lenl){
		df=rbind(df,l[[i]])
	}
	rownames(df)=1:nrow(df)
	return(df)
}

## get x closest peaks to a gene that is in the same TAD as the peak
getClosest <- function(mat,gene,topx=1,minTSS=0,maxDist=1e100,noProm=0){ 
	
	e1=gene
	e2="."
	e3=0
	e4=0
	e5="."
	ef=data.frame(e1,e2,e3,e4,e5)
	enames=c("gene","cond","id","dist","prom")
	colnames(ef)=enames

	tmp=mat[which(mat[,1] == gene),]



	if(nrow(tmp) == 0){
		return(ef)
	}

	## keep only those within max distance
	geneMat=tmp[tmp[,4] < maxDist,]
	
	## just catch too large values here
	if(topx > nrow(geneMat)){
		topx <- nrow(geneMat)
	}

	## order matrix from smallest to biggest distance
	rom=geneMat[order(geneMat[,4]),]

	## implement minTSS and noTSS/noPromoter
	xs=1
	if(minTSS){
		start=which(rom[,4] > minTSS)
		if(length(start) > 0){
			rom=rom[start,]
		}else{
			return(ef)
		}
	}
	if(noProm){
		keep=which(rom[,5] == ".")
		if(length(keep) > 0){
			rom=rom[keep,]
		}else{
			return(ef)
		}
	}
	if(topx > nrow(rom)){
		topx=nrow(rom)
	}

	names(rom)=enames
	return(rom[1:topx,])
}

## table a has all hoxd13 peaks found in the same TAD as a gene
if(0){
	fulltable=read.table("fulltable")

	## a vector with gene ids to check
	c4ids=read.table("cluster4inTAD.gene_ids")

	aout <- list()
	j <- 0
	for(i in as.vector(c4ids[,1])){
		j <- j+1
		aout[[j]] <- getClosest(fulltable,as.character(i),1)
	}
	## make a dataframe out of it
	df <- listToDataFrame(aout)

}

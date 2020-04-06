#!/usr/bin/Rscript


# generate_figure_5J.R -  R script
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
DPATH="../../Fdata/Figure5J/"
SPATH="../../src/"
source(paste(SPATH,"colors.R",sep=""))
source(paste(SPATH,"helper_functions_Fig5J.R",sep=""))
library(svglite)

## filter on minimum mean raw signal in spdh and wt in peak
rawfilter <- function(thresx=-1,inset='indices_to_check',wt.mean='mean_values_of_cond1',spdh.mean='mean_values_of_cond2'){
	if (thresx == -1){
		return(inset)
	}
	si <-which(spdh.mean[inset] > thresx)
	wi <-which(wt.mean[inset] > thresx)
	kb <- intersect(si,wi)
	return(inset[kb])
}


####################################################
##
## READ in rpmbp data                          
##
####################################################

wt   ='L18072' # wt
spdh1='L18074' # het
spdh2='L18076' # hom

spdh=spdh2
reg=11
nfsel=2
if( spdh == "L18076"){
	reg=22
	nfsel=3
}

## normalizing factors from spikeins
##      wt       het        hom
nf=c(1.000752,0.904761,0.766738)

do.norm=1
do.log=T

## these have the 200nt extension and correspond now to what is seen in the browser
wl<-list()
j <- 1
wl[[j]] <- list()
wl[[j]][[1]] <- read.table(paste(DPATH,'/rpmbp_HOXD13_p10_2_centerPeak_',wt,'_25bins_ext200',sep=''),skip=1)
wl[[j]][[2]] <- read.table(paste(DPATH,'/rpmbp_HOXD13_p10_2_centerPeak_',spdh1,'_25bins_ext200',sep=""),skip=1)
wl[[j]][[3]] <- read.table(paste(DPATH,'/rpmbp_HOXD13_p10_2_centerPeak_',spdh2,'_25bins_ext200',sep=""),skip=1)


## add the all regions peaks
wl[[j]][[4]] <- read.table(paste(DPATH,'/rpmbp_all_regions',reg,'_centerPeak_',wt,'_25bins_ext200',sep=""),skip=1)
wl[[j]][[5]] <- read.table(paste(DPATH,'/rpmbp_all_regions',reg,'_centerPeak_',spdh1,'_25bins_ext200',sep=""),skip=1)
wl[[j]][[6]] <- read.table(paste(DPATH,'/rpmbp_all_regions',reg,'_centerPeak_',spdh2,'_25bins_ext200',sep=""),skip=1)
names(wl[[j]]) <- c("WT","HET","HOM","WT.H3K27ac","SPDH.Het_H3K27ac","SPDH.Hom_H3K27ac")

## now we read in the ones with the window size 1kb round the peaks
j <- j+1

wl[[j]]=list()
## the first 2 are also using reads with 200bp expansion
win=1000
wl[[j]][[1]] <- read.table(paste(DPATH,"/rpmbp_HOXD13_p10_2_centerPeak_window_",win,"_centerPeak_",wt,"_25bins",sep=""),skip=1)
wl[[j]][[2]] <- read.table(paste(DPATH,"/rpmbp_HOXD13_p10_2_centerPeak_window_",win,"_centerPeak_",spdh1,"_25bins",sep=""),skip=1)
wl[[j]][[3]] <- read.table(paste(DPATH,"/rpmbp_HOXD13_p10_2_centerPeak_window_",win,"_centerPeak_",spdh2,"_25bins",sep=""),skip=1)
names(wl[[j]]) = c("WT","HET","HOM")
#############

fulltable.f   <- read.table(paste(DPATH,"/fulltable2",sep=""))
fulltable=fulltable.f[,1:5]
allGenes    <- as.matrix(unique(fulltable[,1]))
c4closest.o <- read.table(paste(DPATH,"/cluster_4_dysregulated.bed_closest_hox_midistTSS_0.e",sep=""))
c4Genes     <- as.matrix(c4closest.o[,2])


dfForTargetGenes <- function(target,fmat=fulltable,topx=1,minTSS=0,maxDist=1e100,noProm=0){
	aout <- list()
    j <- 0
    for(i in as.vector(target[,1])){
        j <- j+1
        aout[[j]] <- getClosest(fmat,as.character(i),topx,minTSS,maxDist,noProm)
    }
    ## make a dataframe out of it
	df=do.call(rbind.data.frame, aout)
	return(df)
}
## get the closest peaks now
allGclosest        <- dfForTargetGenes(allGenes,topx=1,minTSS=0,maxDist=1e12)
allGclosest_noProm <- dfForTargetGenes(allGenes,topx=1,minTSS=1000,maxDist=1e12,noProm=1)

fl=list(allGclosest,allGclosest_noProm)
names(fl)=c("allGclosest","allGclosest_noProm")
par(mar=c(9,5,8,1),cex=1.2,cex.main=0.9,mfrow=c(1,4))
pplot <- function(thres=-1,FUN = wilcox.test){
	outl=list()
	j=0
	on=c()
	#for(i in 1:length(fl))
	for(i in 1){
		v=fl[[i]]
		vl=labels(fl)[i]
		############ get PEAK ids ######################

		## the indices and sets to use
		spdh.d=which(v[,2] == "wtUP")
		spdh.u=which(v[,2] == "spdhUP")
		rest.a=(1:nrow(v))[-c(spdh.d,spdh.u)]

		## get the peak ids only 
		spdh.d.p=as.numeric(unique(v[spdh.d,3]))
		spdh.u.p=as.numeric(unique(v[spdh.u,3]))
		rest.a.p=as.numeric(unique(v[rest.a,3]))

		################################################

		cols=3:27
		win=2 ## corresponds here to 1kb win. No window is win=1
		wt.ctrl =rowMeans(wl[[win]]$WT[,cols])
		spdh.het=rowMeans(wl[[win]]$HET[,cols])
		spdh.hom=rowMeans(wl[[win]]$HOM[,cols])

		wt.h3k27 = rowMeans(wl[[1]]$WT.H3K27ac[,cols])
		het.h3k27 = rowMeans(wl[[1]]$SPDH.Het_H3K27ac[,cols])
		hom.h3k27 = rowMeans(wl[[1]]$SPDH.Hom_H3K27ac[,cols])

		## apply normalization factor ###################
		if(do.norm){
			wt.ctrl=wt.ctrl*nf[1]
			spdh.het=spdh.het*nf[2]
			spdh.hom=spdh.hom*nf[3]

			 wt.h3k27=wt.h3k27*nf[1]
            het.h3k27=het.h3k27*nf[2]
            hom.h3k27=hom.h3k27*nf[3]

		}
		pseudo = 0.00001## for log
		fc.all_het=log2((het.h3k27+pseudo)/(wt.h3k27+pseudo))
		fc.all_hom=log2((hom.h3k27+pseudo)/(wt.h3k27+pseudo))
		ic=1:length(fc.all_het)
		fc.all_het.pf=rawfilter(1,inset=ic,wt.mean=wt.h3k27,spdh.mean=het.h3k27)
		fc.all_hom.pf=rawfilter(1,inset=ic,wt.mean=wt.h3k27,spdh.mean=hom.h3k27)
		

		#################################################

		### filter for rpmbp for het
		spdh.u.pf <- rawfilter(thres,inset=spdh.u.p,wt.mean=wt.ctrl,spdh.mean=spdh.het)
		spdh.d.pf <- rawfilter(thres,inset=spdh.d.p,wt.mean=wt.ctrl,spdh.mean=spdh.het)
		rest.a.pf <- rawfilter(thres,inset=rest.a.p,wt.mean=wt.ctrl,spdh.mean=spdh.het)

		j=j+1  ## j=1
		outl[[j]]=log2((spdh.het[spdh.u.pf]+pseudo)/(wt.ctrl[spdh.u.pf]+pseudo))
		j=j+1  ## j=2
		outl[[j]]=log2((spdh.het[spdh.d.pf]+pseudo)/(wt.ctrl[spdh.d.pf]+pseudo))
		j=j+1
		outl[[j]]=log2((spdh.het[rest.a.pf]+pseudo)/(wt.ctrl[rest.a.pf]+pseudo))

		spdh.u.pf <- rawfilter(thres,inset=spdh.u.p,wt.mean=wt.ctrl,spdh.mean=spdh.hom)
		spdh.d.pf <- rawfilter(thres,inset=spdh.d.p,wt.mean=wt.ctrl,spdh.mean=spdh.hom)
		rest.a.pf <- rawfilter(thres,inset=rest.a.p,wt.mean=wt.ctrl,spdh.mean=spdh.hom)

		j=j+1  ## j=1
		outl[[j]]=log2((spdh.hom[spdh.u.pf]+pseudo)/(wt.ctrl[spdh.u.pf]+pseudo))
		j=j+1  ## j=2
		outl[[j]]=log2((spdh.hom[spdh.d.pf]+pseudo)/(wt.ctrl[spdh.d.pf]+pseudo))
		j=j+1
		outl[[j]]=log2((spdh.hom[rest.a.pf]+pseudo)/(wt.ctrl[rest.a.pf]+pseudo))

		on=c(on,paste(vl,c("up-Het","down-Het","rest-Het","up-Hom","down-Hom","rest-Hom"),sep="_"))
	}
	names(outl)=on


	yl=c(-1.6,1.3) ## ylim

	if(do.norm){
	#	yl=c(-2,2) ## restrict ylim range when using normalization
	}

	test.used=FUN(outl[[5]],outl[[6]])$method
	thet <- signif(FUN(outl[[2]],outl[[3]],alternative="less")[[3]],3)
	thom <- signif(FUN(outl[[5]],outl[[6]],alternative="less")[[3]],3)

	o8=outl[1:3]
	o8[[4]]=fc.all_het[fc.all_het.pf]
	o8[[5]]=fc.all_het
	o8[[6]]=outl[[4]]
	o8[[7]]=outl[[5]]
	o8[[8]]=outl[[6]]
	o8[[9]]=fc.all_hom[fc.all_hom.pf]
	o8[[10]]=fc.all_hom

	vl2=as.vector(unlist(lapply(o8,length)))
	on2=paste(vl2,c("up-Het","down-Het","rest-Het","all_H3K27ac_het_filter","all_H3K27ac_het","up-Hom","down-Hom","rest-Hom","all_H3K27ac_hom_filter","all_H3K27ac_hom"),sep="_")
	
	svglite("Figure5J.svg",8,10)
	par(mar=c(15,4,10,2))
	bpb1<-boxplot(o8,
                  col=c("red3","royalblue2","grey","grey","grey","red3","royalblue2","grey","grey","grey"),
                  names=on2,
                  main=paste("closest HOXD13 peak in TAD\nThreshold :",thres,"\nhet-pv = ",thet,"\nhom-pv = ",thom),
                  las=2,ylab="log2 H3K27ac in HOXD13 peak FC Mutant/Wildtype",
                  bty="n",
                  ylim=yl,
                  outline=F
    )

	dev.off()
	return(o8)
}

## make the plot
ret<-pplot(0.01)

#!/usr/local/bin/Rscript

# labbook_analysis_fuctions.R -  R script
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

library(Rtsne)
library(circlize)
library(FNN)
library(igraph)
library(magrittr) ## allows usage of %>% to chain functions, e.g. a %>% sqrt () is the same as sqrt(a)


mycolors=c("#a65628","#c49a6c","#add8e6","#00008b","#377eb8","#ff7f00","#fbb040","#f9ed32","#4daf4a","#f881bf","#e41a1c")
myextcols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000')
mycolors=myextcols
mycolors2 <- c(mycolors,"darkolivegreen1","darkolivegreen3","darkolivegreen4")
mycolors2[8] <- "black"


## now the DBD comparison set to see what we get 
DBD=paste(DPATH,"/lambert/2018_AddAlignments/",sep="")
##          1         2       3       4
dbdsets=c("Homeobox","bHLH","bZIP","Forkhead")
## selected tf sets from the paper
##          1        2          3         4       5       6            7                8           9              10    11      12     13       14      15             16         17
#tfsets=c("C2H2_ZF","KRAB","Homeodomain","bHLH","bZIP","Forkhead","Nuclear_receptor","HMG_Sox","Homeodomain__POU","GATA","SMAD","Runt","Unknown","TRIM","ARID_BRIGHT","AT_hook","C2H2_ZF__Homeodomain","C2H2_ZF__Myb_SANT","CENPB","CUT__Homeodomain","DM","E2F","Ets","Grainyhead","MADS_box","Myb_SANT","Rel","RFX","SAND","T-box","TBP","THAP_finger")
##                  18             19            20            21  22     23    24              25      26        27     28    29     30     31      32
## all tf sets from the paper
alltfsets=c("AP-2","ARID_BRIGHT","ARID_BRIGHT__RFX","AT_hook","BED_ZF","bHLH","Brinker","bZIP","C2H2_ZF","C2H2_ZF__AT_hook","C2H2_ZF__BED_ZF","C2H2_ZF__Homeodomain","C2H2_ZF__Myb_SANT","CBF_NF-Y","CCCH_ZF","CENPB","CG-1","CSD","CSL","CUT__Homeodomain","CxxC","CxxC__AT_hook","DM","E2F","EBF1","Ets","Ets__AT_hook","FLYWCH","Forkhead","GATA","GCM","Grainyhead","GTF2I-like","HMG_Sox","Homeodomain","Homeodomain__Paired_box","Homeodomain__POU","HSF","IRF","KRAB","MADF","MADS_box","MBD","MBD__AT_hook","MBD__CxxC_ZF","mTERF","Myb_SANT","Ndt80_PhoG","NFX","Nuclear_receptor","p53","Paired_box","Pipsqueak","Prospero","Rel","RFX","Runt","SAND","SMAD","STAT","T-box","TBP","TCR_CxC","TEA","THAP_finger","Unknown")

## more visible 
tfsets=c("C2H2_ZF","KRAB","Homeodomain","bHLH","bZIP","Forkhead","Nuclear_receptor","HMG_Sox","Homeodomain__POU","GATA","SMAD","Runt","Unknown","TRIM","ARID_BRIGHT","AT_hook")
#                           17                  18             19            20         21  22     23    24              25      26        27     28    29     30     31 
tfsets=c(tfsets,"C2H2_ZF__Homeodomain","C2H2_ZF__Myb_SANT","CENPB","CUT__Homeodomain","DM","E2F","Ets","Grainyhead","MADS_box","Myb_SANT","Rel","RFX","SAND","T-box","TBP")
##                 32               33                        34
tfsets=c(tfsets,"THAP_finger","Homeodomain__Paired_box","Paired_box")

## function to determine diAA aboundance
diA <- function(x){
    res=c()
    for(i in 1:19){
        for(j in 2:20){ 
            if(i < j){
                res=c(res,x[i]+x[j])
            }
        }   
    }
    return(res)
}


## function to read sets of files into a list processed from before

setP=paste(DPATH,"/Figure7/protparam/",sep="")
readAAm <- function(sel=1,insets=tfsets,P=setP,addf='_protparam_matrix_IDR_new_filtered',filtcol=c(5:17,19:23,25:31),AAcol=2:21){

	mcolsindex=c()
	mrowsnames=c()
	inputl <-list()
	inputl2<-list()
	inl=0

	for(nset in insets[sel]){
		inl = inl+1
		fin=paste(P,nset,addf,".csv",sep="")

		## read in matrix and header
		a=read.table(fin)
		header=scan(fin,nlines=1,what="character")
		colnames(a)=header

		## need to exclude PYL and Sel since they are 0, also excluding the non numerical values. Col32 is the actual sequence inspected
		if(filtcol[1] == 0){
			b=a
		}else{
			b=as.matrix(a[,filtcol])
		}
		rownames(b)=a[,1]
		#bs=scale(b)
		## add new set to colindex
		mcolsindex=c(mcolsindex,rep(inl,dim(b)[1]))
		mrowsnames=c(mrowsnames,as.vector(a[,1]))

		## maybe scaling should be done after all sets are read in?
		## we need to transpose here since the unlisting of a matrix is by column.
		inputl[[inl]]=t(b)
		inputl2[[inl]]=t(a)
		#pp(bs,nset)
	}
	outm=matrix(unlist(inputl),ncol=dim(inputl[[1]])[1],byrow=T)
	fullm=matrix(unlist(inputl2),ncol=dim(inputl2[[1]])[1],byrow=T)
	
	colnames(outm)=c(rownames(inputl[[1]]))
	rownames(outm)=mrowsnames
	colnames(fullm)=c(rownames(inputl2[[1]]))

	outl=list(mcolsindex,mrowsnames,outm,fullm,insets[sel])

	if(AAcol[1] == 0){
		outl[[6]]=1
	}else{
	## get diAA frequencies
		header1<-colnames(b)[AAcol]
		la=length(AAcol)
		getdin<-c()
		for(i in 1:(la-1)){
			for(j in 2:la){ 
				if(i < j){
					getdin=c(getdin,paste(header1[i],header1[j],sep=""))
				}
			}
		}
		outl[[6]] <- t(apply(outl[[3]][,AAcol],1,diA))
		colnames(outl[[6]]) <- getdin
	}
	names(outl)=c("groupindex","genes","matrixvalues","fullmatrix","sets","diAA")
	return(outl)
}

## BIC says optimal cluster number is 8 which is in accordance of AIC as well
bicplot <- function(pcao,plotme=T,dims.use=10,biconly=F,seed=999){
	clist<-list();
	j=0
	for(nclust in 2:22){
		j=j+1;
		set.seed(seed)
		km<-kmeans(pcao$x[,1:dims.use],center=nclust);
		clist[[j]]=BIC2(km)
	}

	## these are the cluster numbers checked
	names(clist)=2:22
	#names(clist)[which(unlist(clist)[seq(2,63,3)] %in% min(unlist(clist)[seq(2,63,3)]))]

	minB=order(unlist(clist)[seq(2,63,3)])[1]+1 ## since we start counting at 2
	if(plotme){

		## first plot bic
		miny = min(unlist(clist))
		maxy =max (unlist(clist))
		plot(2:22,unlist(clist)[seq(2,63,3)],"l",xlab="Number of clusters",ylab="[AU]",main=paste("Minimum BIC =",minB),xlim=c(2,22),col="cornflowerblue",ylim=c(miny,maxy))
		if(biconly == F){
			points(2:22,unlist(clist)[seq(3,63,3)],col="grey",type="b",pch=17)
			lines(2:22,unlist(clist)[seq(1,63,3)])
			legend("topright",c("AIC","BIC","AICc"),col=c("black","cornflowerblue","grey"),lwd=c(2,2,2),pch=c(-1,-1,17))
		}else{
			legend("topright",c("BIC"),col=c("cornflowerblue"),lwd=c(2))

		}
	}
	return(clist)
}

## the fit function is for testing and essential
bicplot.fit <- function(pca,s=546,offset=1,pm=T,dims.u=11){ ## the used dims are selected in a way thay we keep more then 80% of the original variance
	clist=bicplot(pca,dims.use=dims.u,biconly=T,seed=s,plotme=pm)
	y=unlist(clist)[seq(2,63,3)]
	x=2:22
	##fit <- lm(log(y)~(log(x)+log(x^2)+log(x^3)))
	fit <- lm(log(y)~log(x))
	yfit <- x^fit$coefficients[2]*exp(fit$coefficients[1])
	##plot(x,y,col="cornflowerblue","l")

	of=offset
	dd <- (abs(y-yfit))[(1+of):(length(x)-of)]
	maxi<-which(dd==max(dd))+of+1
	if(pm){
		lines(x,yfit,col="red")
		points(x,yfit,col="red",pch=16)
		abline(v=maxi,lty=2,col="grey")
		abline(h=y[(maxi)],lty=2,col="grey")

		#########################################################
		lines(c(maxi,maxi),c(y[(maxi-1)],yfit[(maxi-1)]),lwd=3)
		text(maxi,yfit[(maxi-1)],paste("k=",maxi),pos=4)
	}
	return(maxi)
}


bicplot.old <- function(pcao,plotme=T,dims.use=5){
clist<-list();
j=0
for(nclust in 2:22){
    j=j+1;
    km<-kmeans(pcao$x[,1:dims.use],center=nclust);
    clist[[j]]=BIC2(km)
}
## these are the cluster numbers checked
names(clist)=2:22
names(clist)[which(unlist(clist)[seq(2,63,3)] %in% min(unlist(clist)[seq(2,63,3)]))]

minB=order(unlist(clist)[seq(2,63,3)])[1]+1 ## since we start counting at 2
if(plotme){
plot(2:22,unlist(clist)[seq(1,63,3)],"l",xlab="Number of clusters",ylab="[AU]",main=paste("Minimum BIC =",minB),xlim=c(2,22))
points(2:22,unlist(clist)[seq(3,63,3)],col="grey",type="b",pch=17)
lines(2:22,unlist(clist)[seq(2,63,3)],col="cornflowerblue")
legend("topright",c("AIC","BIC","AICc"),col=c("black","cornflowerblue","grey"),lwd=c(2,2,2),pch=c(-1,-1,17))
}
return(clist)
}

## smooth it by repetition
bicrep <- function(dend,du=10,reps=1000){
	clist<-bicplot(dend,plotme=F,dims.use=du)
	final=list()
	for(j in 1:3){
		final[[j]]=unlist(clist)[seq(j,63,3)]
	}

	for(i in 2:reps){
		clist<-bicplot(dend,plotme=F,dims.use=du)
		for(j in 1:3){
			final[[j]]=final[[j]]+unlist(clist)[seq(j,63,3)]
		}
	}
	for(j in 1:3){
		final[[j]]=final[[j]]/reps
	}
	return(final)
}


## function to process the input matrices and filter on AA di or single frequencies		
process.mat <- function(inset=tf.list,type="sAA",thresAA=0.2){
	## outmat
	om <- data.frame()

	if(type=="dAA"){
		om <- data.frame(which(inset[[6]] > thresAA,arr.ind=T))
		om[,2]<-getdin[om[,2]]
	}else{
		## outmat
		om <- data.frame(which(inset[[3]][,2:21] > thresAA,arr.ind=T))
		om[,2]<- colnames(inset[[3]][,2:21])[om[,2]]
	}
	## which of those are left now and we can use for clustering
	#setnew<-unique(om[om[,1]> 76,1])
	## om can have the same gene in twice if more than one AA has fraction > threshold
	setnew<-sort(unique(om[,1])) ## we sort so we have the same row order as in the original matrices

	## now make the initial approach again and run pca and clustering on reduced set and full set including diAA frequencies
	nmat=inset[[3]][setnew,]
	#nmat.e=cbind(lout[[3]],diAmat)
	#nmat.e.red <- nmat.e[setnew,]

	dups=which(duplicated(om[,1])==T)
	## dups has the indices of duplicated rowids
	omtmp=om
	
	if(length(dups) > 0){
	for(i in dups){
		## now we get all row ids which have the same id as in dups and then fuse the results to a new singleton id
		tk<-which(om[,1]%in%om[i,1])
		omtmp[tk[1],2]=paste(om[tk,2],collapse=",")
	}

	ot2 <- omtmp[-dups,]
	}else{
		ot2 <- omtmp
	}
	ot3<- ot2[order(ot2[,1]),]

	tmat=nmat
	if(!is.matrix(tmat)){
		tmat=matrix(tmat,ncol=length(tmat),nrow=1)
		rownames(tmat)=rownames(om)
		colnames(tmat)=colnames(inset[[3]])
	}

	## ok up to here it works nicely 
	rown2=paste(addbraces(ot3[,2]),rownames(tmat))
	params=c(type,thresAA,dim(inset[[3]])[1],length(rown2))
	names(params)=c("filter","thres","#in","#out")
	outl=list(tmat, rownames(tmat),ot3[,2],rown2,setnew,params)
	names(outl)=c("matrix","geneids","AA","fused","indices","parameters")

	return(outl)
}

getdend <- function(tmat=dmat,dims.use=10,nclust=5,md="complete",scree=F,seed=999){

## scale and center
test.mat=scale(tmat)
## PCA
set.seed(seed)
pca<-prcomp(test.mat,center=F,scale.=F)

## making scree plot
if(scree){
	myscree(pca,pca=T,cs=T)
}

## how many of the 25 PCs to use
#dims.use=10 ## cover ~80% of variance from the data
#dims.use=5 ## cover ~50% of variance from the data


#nclust=5
set.seed(seed)
km<-kmeans(pca$x[,1:dims.use],center=nclust)

t1=pca$x[,1:dims.use]
mytext=unlist(strsplit(rownames(t1)," "))[seq(1,2*dim(t1)[1],2)]
rownames(t1) = 1:dim(t1)[1]

## dist <- calculates distances between rows of a matrix, default is euclidian
## hclust makes hierachical clustering based on a clustering method. default is 'complete' linkage! we can also try ward and UPGMA and single linkage clustering
## adjust method here. Either ward.D/D2, single, complete, average
dend <- as.dendrogram(hclust(dist(t1),method=md))

## getting the rowIDs as numeric values as they go from left to right in the dendrogram!!!
dorder <- as.numeric(labels(dend))

## getting the colors from kmeans
mykmeanscol=as.numeric(km$cluster[dorder])

## setting the color of each leaf now to the kmeans cluster number assigned!!!
labels_colors(dend)=mycolors2[mykmeanscol]
labels(dend)=rownames(pca$x)[dorder]             ## this should not be necessary
ol <- list(dend,dorder,pca,km,c(dims.use,nclust))
labels(ol[[5]])=c("pcs.used","k-clusters")
names(ol)=c("dendrogram","dorder","pca","km","parameters")
return(ol)
}


## this function reorders the top branches based on their cluster proximity, paper version only !
hkm <- function(dall=idrs_dall,dims.use=10,nclust=7,seed=999,both=F,rorder=1){
    dlist <-list()
    for(i in 1:nclust){
        indic=which(dall$km$cluster ==i)
        matin=dall$pca$x[,1:dims.use]

        ## we keep as rownames the orignal matrix numbers and just add them later again
        rownames(matin)=1:(dim(matin)[1])
        set.seed(seed)
        dlist[[i]]=as.dendrogram(hclust(dist(matin[indic,])))
        ## color leaves and rename
        #labels_colors(dlist[[i]])=mycolors[i]
        labels_colors(dlist[[i]])=i  ## assign just the cluster number as the color so we can later just choose a color we like
        #   dorder <- as.numeric(labels(dlist[[i]]))
        #   labels(dlist[[i]])=rownames(matin)[dorder]
    }

	if(length(rorder) > 1){   ## this should become the default!!!! anyhow we keep it for consitency until paper is out
		## we need to automate this but for now it shall be ok
		if(nclust == 7){ ## we need to call this one if having rearranged with remat and we want to get the same dendrogram as in the paper
			## this is of course for comparison reasons only but nice to have
			## observe that the numbers in the assembly here correspond to the relabeled dendrogram numbers from the paper version
			## 
			df1 <- merge(dlist[[1]],dlist[[7]])
			df1 <- merge(df1,dlist[[6]])
			df1 <- merge(df1,dlist[[5]])
			df1 <- merge(df1,merge(dlist[[4]],dlist[[3]]))   ## cluster 4 new was cluster2 orig, 
			df1 <- merge(df1,dlist[[2]])                     ## cluster 2 new was cluster6 orig
			## check the values below after both=F to see that this is indeed true
		}else{
			ldlist=length(dlist)
			df1=merge(dlist[[1]],dlist[[nclust]])
			for(i in c((nclust-1):2)){
				df1 <- merge(df1,dlist[[i]])
			}
		}
	}else{## to get in here reorder is a vector of length one. Original paper version
		## merge all of them
		if(nclust == 4){
			df1<- merge(merge(dlist[[1]],dlist[[4]]),merge(dlist[[2]],dlist[[3]]))
		}else if(nclust == 7){ ### this is the original dendrogram
			if(both==F){
				df1<-merge(merge(merge(merge(merge(dlist[[1]],dlist[[7]]),dlist[[5]]),dlist[[4]]),merge(dlist[[2]],dlist[[3]])),dlist[[6]])
			}else{
				df1<-merge(merge(merge(merge(merge(dlist[[1]],dlist[[7]]),dlist[[5]]),dlist[[4]]),merge(dlist[[3]],dlist[[2]])),dlist[[6]])
				}
		}else{
			seen=c(1)
			df1<-dlist[[1]]
			for(j in 2:nclust){
				df1<-merge(df1,dlist[[j]])
			}
		}
	}
    return(df1)
}


#!/usr/local/bin/Rscript

# circplot_fun_re.R_new -  R script
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


mtrack <-function(mycol=ccol,ymax=0.2,yh=0.1,th=0.07,mar=c(0)){
	if(length(mar) > 1){
		circos.par(track.margin=mar)
	}
	n=length(mycol);
	circos.track(ylim = c(0, ymax), panel.fun = function(x, y) {
					 circos.rect(0:(n-1),0,1:n, yh,col=mycol,border=NA)
},bg.border=NA,track.height=th)
}


hlgene <-function(gene="HOXD",pm=pmat,od=df1_labels){
	mytcol=rep(1,length(pm[[2]]))
	mtrack(mycol=mytcol[od])
	mytcol[grep(gene,pm[[2]])]=3
	mtrack(mycol=mytcol[od])
}

## good
myAAcols=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

myAA=sort(c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"))

myAAcols[14]="skyblue1"

## a matrix that assigns to each AA one color !!!
AAcolF=cbind(myAA,myAAcols[1:20])
AAcolF[1,2] = colors7[1]
AAcolF[14,2] = colors7[6]

circplot <-function(pmat=pmat.tf2_19,tfl=tf.list,pos_set=pos_set10,nfam=8,nclust=7,dims.use=10,plotme=T,seed=999,rorder=1,retype='NA',remat=matrix()){
	## select matrix we want to use
	#pmat=pmat.tf2_19
	#names(pmat)
	#    1         2       3      4        5         6
	# "matrix"  "geneids" "AA" "fused" "indices" "parameters"

	## orignal input matrix for pca 
	dmat=pmat[[1]] 
	## get_dendrogram
	idrs_dall <-getdend(dmat,dims.use=dims.use,nclust=nclust,md="ward.D2",seed=seed)

	## I should reassign cluster ids here already. Will make all things easier
	gene_ids_cluster=as.numeric(idrs_dall[[4]]$cluster)
	gene_ids_cluster.old=gene_ids_cluster

	## this is a new routine which should make it more easy to rearrange cluster ids
	## it also resets the original kmeans ids 
	## reorder by mapping matrix given to function
	if(nrow(remat) > 1){
		onames=names(idrs_dall[[4]]$cluster)
		gene_ids_cluster <- rep(1,length(gene_ids_cluster.old))
		## reordering is done by checking which of the original paper clusters are most caught in the new ids. Then the new ids are remapped to the old ones!
		## check lines 61 to 102 in labbook_takeout_ZF.R_new
		for(i in 1:nrow(remat)){	                      ## current  -> New
			gene_ids_cluster[which(gene_ids_cluster.old == remat[i,1])]= remat[i,2]
		}
		## reset now the values
		idrs_dall[[4]]$cluster <- gene_ids_cluster
		names(idrs_dall[[4]]$cluster) <- onames
	}

	if(retype == 'ZF' & nclust == 7){
		onames=names(idrs_dall[[4]]$cluster)
		gene_ids_cluster <- rep(1,length(gene_ids_cluster.old))
		## reordering is done by checking which of the original paper clusters are most caught in the new ids. Then the new ids are remapped to the old ones!
		## check lines 61 to 102 in labbook_takeout_ZF.R_new

		gene_ids_cluster[which(gene_ids_cluster.old == 2)]= 1
		gene_ids_cluster[which(gene_ids_cluster.old == 7)]= 2
		gene_ids_cluster[which(gene_ids_cluster.old == 5)]= 3
		gene_ids_cluster[which(gene_ids_cluster.old == 3)]= 4
		gene_ids_cluster[which(gene_ids_cluster.old == 1)]= 5
		gene_ids_cluster[which(gene_ids_cluster.old == 6)]= 6   ## stays the same
		gene_ids_cluster[which(gene_ids_cluster.old == 4)]= 7

		## reset now the values
		idrs_dall[[4]]$cluster <- gene_ids_cluster
		names(idrs_dall[[4]]$cluster) <- onames
	}
	if(retype == 'KRAB' & nclust == 6){
        onames=names(idrs_dall[[4]]$cluster)
        gene_ids_cluster <- rep(1,length(gene_ids_cluster.old))
        ## reordering is done by checking which of the original paper clusters are most caught in the new ids. Then the new ids are remapped to the old ones!
        ## check lines 61 to 102 in labbook_takeout_ZF.R_new
                                                     #new   old
        gene_ids_cluster[which(gene_ids_cluster.old == 1)]= 5 ## 100% from what is left
        gene_ids_cluster[which(gene_ids_cluster.old == 2)]= 4 ## 95%
        gene_ids_cluster[which(gene_ids_cluster.old == 3)]= 6 ## 98 % 
        gene_ids_cluster[which(gene_ids_cluster.old == 4)]= 1 ## 98% recovery
		#gene_ids_cluster[which(gene_ids_cluster.old == 5)]= 7  ## 61 % recovery 
        gene_ids_cluster[which(gene_ids_cluster.old == 5)]= 3  ## 61 % recovery  ## since there are only 6 clusters we use cluster 3 now 
		                                                       ## as old cluster 7 recovery!!!
        gene_ids_cluster[which(gene_ids_cluster.old == 6)]= 2   ## 93% recover 

		## we lost old cluster 3 here due to having one cluster less and the fact that lold cluster 3 is also distributed to 3 different new 
		##clusters with 16,26 and 47 to new clusters 2,4,5

        ## reset now the values
        idrs_dall[[4]]$cluster <- gene_ids_cluster
        names(idrs_dall[[4]]$cluster) <- onames
    }
#	if(retype == 'ORIG' & nclust == 7){
#		onames=names(idrs_dall[[4]]$cluster)
#		gene_ids_cluster <- rep(1,length(gene_ids_cluster.old))
#		## reordering is done by checking which of the original paper clusters are most caught in the new ids. Then the new ids are remapped to the old ones!
#		## check lines 61 to 102 in labbook_takeout_ZF.R_new
#
#		gene_ids_cluster[which(gene_ids_cluster.old == 2)]= 1
#		gene_ids_cluster[which(gene_ids_cluster.old == 7)]= 2
#		gene_ids_cluster[which(gene_ids_cluster.old == 5)]= 3
#		gene_ids_cluster[which(gene_ids_cluster.old == 3)]= 4
#		gene_ids_cluster[which(gene_ids_cluster.old == 1)]= 5
#		gene_ids_cluster[which(gene_ids_cluster.old == 6)]= 6   ## stays the same
#		gene_ids_cluster[which(gene_ids_cluster.old == 4)]= 7
#
#		## reset now the values
#		idrs_dall[[4]]$cluster <- gene_ids_cluster
#		names(idrs_dall[[4]]$cluster) <- onames
#	}
#	return(idrs_dall[[4]]$cluster)



	#names(idrs_dall)
	#    1           2           3           4                5
	#"dendrogram" "dorder"     "pca"        "km"         "parameters"

	## first has the rowids in proper order 
	#tmp1<-idrs_dall[[3]]$x
	#rownames(tmp1)=1:dim(tmp1)[1]
	## get dendrogram
	#dend <- as.dendrogram(hclust(dist(tmp1),method="ward.D2"))

	## second uses the gene ids assigned
	#tmp2=idrs_dall[[3]]$x
	## get dendrogram
	#dend2 <- as.dendrogram(hclust(dist(tmp2),method="ward.D2"))

	## now we compare if the labels are the same!
	## leafs        orig rownames   numeric ids
	#labels(dend2)==rownames(dmat)[as.numeric(labels(dend))]

	### matches so this is correct!

	## lets check now if the cluster ids are correct!
	#names(idrs_dall[[4]]$cluster)==rownames(tmp2)  ## that fits

	## thus => 
	## each gene in orig order gets its cluster number assigned

###################################
###################################  Make the inner dendrogram and arrange the clusters in the order we want
###################################  This determines the final output appearance 
###################################  The cluster IDs should already be reassigned optimally
###################################  Check labbook_analysis_functions.R and there the hkm method
###################################  The cluster ids from kmeans are in object idrs_dall$km$cluster
###################################  Check labbook_takeout_ZF.R_new for part -> old new id comparison to see how to reassign the ids for the ZF analysis

	## this just takes the kmeans clusters, makes a tree for each cluster
	## and then puts the trees togther, the order can be refined within the function -> th
	df1=hkm(idrs_dall,nclust=nclust,dims.use=dims.use,rorder=rorder)

###################################
###################################  Make the inner dendrogram and arrange the clusters in the order we want
###################################  
###################################

	## make a new tree here that fits with our image for the paper

	#plot(df1) ## shows the new dendrogram 
	#labels(df1) ## need to be the ids of the genes per cluster

	## the genes in cluster 1 must be in the very front etc ...
	#which(gene_ids_cluster == 1)  ## red color => fits all correct
	#which(gene_ids_cluster == 2)  ## green color => fits all correct
	#which(gene_ids_cluster == 5)  ## orange color => fits all correct

	## the labels from left to right in the dendrogram are
	df1_labels <- as.numeric(labels(df1))
	df1_gene_clusters <- gene_ids_cluster[df1_labels]  ## fits
	## clusters are merged as
	# c(1, 7, 5, 4, 2, 3, 6)
	cs=c() ## cluster sizes
	for(i in c(1, 7, 5, 4, 2, 3, 6)){
		cs=c(cs,length(which(sort(gene_ids_cluster) == i)))
	}
	#
	### ids in each cluster in dendrogram in linear order are the same as the ids found in the
	### kmeans clustering => new tree with leafs has correct ids per cluster assigned !!!
	#j=0;
	#ts=1
	#for(i in c(1, 7, 5, 4, 2, 3, 6)){
	#	j=j+1
	#	tt <- sort(which(gene_ids_cluster == i))
	#	from=ts
	#	to=from+cs[j]-1
	#	ts=ts+cs[j]
	#
	#print(sort(df1_labels[(from:to)])==tt )
	#}
	#
	###### we can take the leaf ids now / which are the indices from the input matrix
	###### and get all other values !
	nsize=length(df1_labels) ## just the length
	## set df1 leaf labels now to gene names
	#labels(df1)=rownames(dmat)[df1_labels]

	## get the hox genes in this list!!!

	#grep(" HOX",labels(df1))

	## get TF family for each idr
	#nfam=8  ## number of families we want to include. eg 7 named ones +1 for the rest to fuse

	famind.idr <-  tfl$groupindex[pmat$indices]
	#famind.dbd <- dbd.list$groupindex[pmat$indices]
	famind.idr1 = famind.idr
	## set lowlevel fams all to 8
	famind.idr1[famind.idr1>=nfam]=nfam


	## and now the ids which are in the positive set of TFs
	posids <- which(pmat[[2]] %in% pos_set[,1])
	revposidtable<-pos_set[which(pos_set[,1] %in%pmat[[2]]),] ## this table has the indices of the full matrix from the positive set
	## some dont have a proper repeat of length 10 and should be taken out !!!!

	## pmat[[2]][posids]
	mytcol=rep(1,nsize)
	mytcol2=rep(1,nsize)
	#mytcol2=rep(1,n)
	#mytcol[posids]=2 ## just a different color

	## get the repeats from IDRs we collected
	repeats <- tf.list$fullmatrix[pmat$indices,c(1,33:34)]

	valid_idrs <- which(as.numeric(repeats[,3])>=10) ## the idrs which have a valid repeat => len >=10
	posids_valid <- intersect(valid_idrs,posids)    ## all idrs with a repeat > 10 and being in the positive set from the data table

	ids_valid_list <-list(posids,posids_valid,valid_idrs)
	names(ids_valid_list) <- c("posids","posids_valid","valid_idrs")

	##
	#posids2 <- posids[which(posids %in%valid)]
	#revposidtable2<-pos_set[which(pos_set[,1] %in%pmat[[2]][posids2]),]
	revposidtable2<-pos_set[which(pos_set[,1] %in%pmat[[2]][posids_valid]),]  ## this is our gold standard table 
	## positve set ids and idr length > 10  !!! 

	## we just for now need to take the posids and the repeat that is annotated 


	## I need to check if the Repeat is actually the same between the published table and the IDR I have so far

	cin=c()
	cout=c()
	for(i in 1:dim(repeats[posids_valid,])[1]){
		ind<-which(revposidtable2[,1] %in% repeats[posids_valid,][i,1])
		if(repeats[posids_valid,][i,2] == revposidtable2[ind,6]){
			cin=c(cin,ind)
		}else{
			cout=c(cout,ind)
		}
	}

	## 41 are ok
	## 3 are different
	## naja

	## we base now the repeats on the table only. Later we can take also what we have in the IDRs

	revposidtable_s <- revposidtable
	reod=c()
	for(i in pmat[[2]][posids]){
		reod=c(reod,which(revposidtable[,1] == i))
	}
	revposidtable_s <- cbind(revposidtable[reod,],posids,repeats[posids,2:3])
	revposidtable_s <- cbind(revposidtable_s,as.vector(revposidtable_s[,6]) == as.vector(revposidtable_s[,8]),as.numeric(revposidtable_s[,9]) >=10)
	colnames(revposidtable_s)[10:11]=c("sameAA","validIDRAA")
	#rownames(revposidtable_s)=1:66
	rownames(revposidtable_s)=1:nrow(revposidtable_s)

	## finally it is correct
	replist=list();j=0                                                                ## j=  1 2 3   
	for(nn in names(table(revposidtable_s[,6]))){ ## for the repeats in the poslist, nn iselem [A,C,S,...]
		j=j+1;
		replist[[j]]=revposidtable_s[which(revposidtable_s[,6]==nn),7]
	}
	names(replist)=names(table(revposidtable_s[,6]))

	## only the valid ones and the AA from the annotated TF
	replist2=list();j=0                                                                ## j=  1 2 3   
	for(nn in names(table(revposidtable_s[,6]))){ ## for the repeats in the poslist, nn iselem [A,C,S,...]
		j=j+1;replist2[[j]]=revposidtable_s[ which(revposidtable_s[,6]==nn & revposidtable_s[,11] == T),7]
	}
	names(replist2)=names(table(revposidtable_s[,6]))


	## only the valid ones and the AA from the IDR we caught
	onlyv <- which(revposidtable_s[,11] == T)
	replist3=list();j=0                                                                ## j=  1 2 3   
	for(nn in names(table(revposidtable_s[onlyv,8]))){ ## for the repeats in the poslist, nn iselem [A,C,S,...]
		j=j+1;replist3[[j]]=revposidtable_s[ which(revposidtable_s[,8]==nn & revposidtable_s[,11] == T),7]
	}
	names(replist3)=names(table(revposidtable_s[onlyv,8]))


	## now lets make the circle plot colors for the IDRS etc
	mytcol=rep("NA",nsize)
	mytcol2=rep("NA",nsize)
	mytcol3=rep("NA",nsize)
	mytcol4=rep("NA",nsize)  ## colors for repeats we called in our IDRs

	## all I have in my data
	#AAcol=cbind(1:length(names(replist)),names(replist))
	#AAcolF ## has a color foreach of the 20 AA


	allAAin=c()
	namesreplist=names(replist)
	for(j in 1:length(replist)){
		cj=which(AAcolF[,1] %in%namesreplist[j])
		mytcol[replist[[j]]]= AAcolF[cj,2]
	}
	allAAin=c(allAAin,namesreplist)
	namesreplist=names(replist2)
	for(j in 1:length(replist2)){
		cj=which(AAcolF[,1] %in%namesreplist[j])
		mytcol2[replist2[[j]]]= AAcolF[cj,2]
	}
	allAAin=c(allAAin,namesreplist)
	namesreplist=names(replist3)
	for(j in 1:length(replist3)){
		cj=which(AAcolF[,1] %in%namesreplist[j])
		mytcol3[replist3[[j]]]= AAcolF[cj,2]
	}
	allAAin=c(allAAin,namesreplist)
	allAAin=sort(unique(allAAin))


	############# COLORS FOR REPEATS ###########
	## get colors for our repeats called in IDRs
	for(j in valid_idrs){
		cj<-which(AAcolF[,1] %in% repeats[j,2])
		mytcol4[j]=AAcolF[cj,2]
	}
	##
	############################################


	## now make are circos plot again
	# df1_labels has the numeric ids in as the leafs in the dendrogram, which are the original indices in pmat

	## this is all for plotting now

	## init plot region
	if(plotme){
	circos.clear()
	circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0),gap.degree=0)
	circos.initialize("foo", xlim = c(0, nsize))

	## depending on the color we supply we get a corresponding circle from outside to inside
	## the colors need to have the same order as the value in the dendrogram is!

	## plot TF family has simply base R colors 1 to 8 for the first 8 families
	ccol=famind.idr1[df1_labels]
	mtrack(mycol=ccol)

	## plot positive set -> colors from mycolors[1:j] where j is the number of different AA repeats found in the data => so A maybe 1, G maybe 2 ... 
	#ccol=mytcol[df1_labels]
	#mtrack(mycol=ccol)

	## plot positive set IDRs with paper TF AA
	#ccol=mytcol2[df1_labels]
	#mtrack(mycol=ccol)

	## plot positive set IDRs with IDR AA found -> 49 are the same, 4 are different
	#ccol=mytcol3[df1_labels]
	#mtrack(mycol=ccol)

	## plot all IRDs with an AA repeat found >= 10
	ccol=mytcol4[df1_labels]
	mtrack(mycol=ccol)

	## colors from hkm function using mycolors[1:7] for the the 7 clusters
	ccol=labels_col(df1)
	mtrack(mycol=ccol)


	## legend for clusters -> the colors come from hkm function !!!!
	#legend("topright",paste(paste("cluster",1:nclust),  addbraces(table(gene_ids_cluster))),pch=16,col=mycolors[1:nclust])
	legend("topright",paste(paste("cluster",1:nclust),  addbraces(table(gene_ids_cluster))),pch=16,col=1:nclust)
	## legend for tf famlies

	legend("topleft",paste(    c(tfl$set[1:(nfam-1)],   paste("REST_",32-nfam+1,"_families",sep="")  ), addbraces(table(famind.idr1))   ),pch=15,col=1:nfam)

	## legend for repeat IDRs
	#legend("bottomright",paste(names(replist),col=mycolors[1:length(replist)],lwd=2)

	### merge AA col from postivIDs with posidtable and the values there 
	#posidsAA<-table(revposidtable[,6])
	#posidsAA2<-table(revposidtable2[,6])
	#AAcol2=cbind(AAcol,rep(0,dim(AAcol)[1]),rep(0,dim(AAcol)[1]))
	#for(i in 1:length(AAcol[,2])){
	#	ind=which(names(posidsAA)%in% AAcol[i,2])
	#	if(length(ind)>0){
	#		AAcol2[i,3]=posidsAA[ind]
	#	}
	#	ind=which(names(posidsAA2)%in% AAcol[i,2])
	#    if(length(ind)>0){
	#        AAcol2[i,4]=posidsAA2[ind]
	#    }
	#}
	#
	## now we have a table finally 

	#allAAin ## sorted AAs from the lists

	to=c(); ## should hold then increasing indices anyways
	for(i in allAAin){ 
		to=c(to,which(AAcolF[,1] %in% i))    
	}


	## make the legend matrix
	nmat=matrix(nrow=length(allAAin),ncol=3)
	rownames(nmat)=allAAin
	for(i in 1:length(allAAin)){
		nn=allAAin[i]
		ind <- which(names(replist)==nn)
		if(length(ind) > 0){
			nmat[i,1]=length(replist[[ind]])
		}else{
			nmat[i,1]=0
		}
	}

	for(i in 1:length(allAAin)){
		nn=allAAin[i]
		ind <- which(names(replist2)==nn)
		if(length(ind) > 0){
			nmat[i,2]=length(replist2[[ind]])
		}else{
			nmat[i,2]=0
		}
	}


	for(i in 1:length(allAAin)){
		nn=allAAin[i]
		ind <- which(names(replist3)==nn)
		if(length(ind) > 0){
			nmat[i,3]=length(replist3[[ind]])
		}else{
			nmat[i,3]=0
		}
	}

	fnl=c()
	for(i in 1:dim(nmat)[1]){
		fnl=c(fnl,paste(rownames(nmat)[i],paste0(nmat[i,],collapse=" ")))
	}

	#legend("bottomleft",fnl,col=AAcolF[to,2],lwd=2)
	#table(repeats[valid_idrs,2])
	rti<-table(repeats[valid_idrs,2])
	fnl=paste(names(rti),as.numeric(rti))
	fnlc=c()
	for(j in names(rti)){
		cj=which(AAcolF[,1] %in% j)
		fnlc=c(fnlc,AAcolF[cj,2])
	}
	legend("bottomright",fnl,col=fnlc,lwd=2)
	}

	inputargs=list(pmat,tfl,pos_set,nfam,nclust,dims.use)
	names(inputargs)=c("pmat","tfl","pos_set","nfam","nclust","dims.use")
	#4 cluster    #5 family,  ## revpostab  ## result of getded
	outlist=list(df1_labels,revposidtable,revposidtable2,gene_ids_cluster,famind.idr,revposidtable_s,idrs_dall,inputargs,ids_valid_list,df1)
	names(outlist)=c("df1_labels","revposidtable","revposidtable2","gene_ids_cluster","famind.idr","revposidtable_s","idrs_dall","inputargs","ids_valid_list","df1")
	return(outlist)
}

## replot the circle from tla file created by the circplot function itself
## in this function we cosmetically just reassign the cluster numbers so they match
## what is found in the paper figure numbering

replot_circle <- function(obj=tla,ppset=4,plotdend=F,cexl=1.2,colsClust=1:7,colsFam=1:8,hlf=F,reassign=1,reordertree=F,hlfg=c(),tffam=T){
	idrs_all <- obj[[7]]
	df1_labels<- obj[[1]]
	nsize=length(df1_labels)
	gene_ids_cluster.old<- obj[[4]]
#	nclust <- krab$idrs_dall$parameters[[2]]
#	dims.use <- krab$idrs_dall$parameters[[1]]

	if(reassign == 1){
		## change cluster ids here but not the dendrogram order!
		gene_ids_cluster <- rep(1,length(gene_ids_cluster.old))

		## cluster   2 -> 4
		## cluster   4 -> 5
		## cluster   5 -> 6
		## cluster   6 -> 2
		gene_ids_cluster[which(gene_ids_cluster.old == 2)]= 4
		gene_ids_cluster[which(gene_ids_cluster.old == 4)]= 5
		gene_ids_cluster[which(gene_ids_cluster.old == 5)]= 6
		gene_ids_cluster[which(gene_ids_cluster.old == 6)]= 2
	
		## those stay the same	
		gene_ids_cluster[which(gene_ids_cluster.old == 1)]= 1
		gene_ids_cluster[which(gene_ids_cluster.old == 7)]= 7
		gene_ids_cluster[which(gene_ids_cluster.old == 3)]= 3
	}else{
		gene_ids_cluster <- gene_ids_cluster.old
	}

	famind.idr <- obj[[5]]
    famind.idr1 = famind.idr

	nfam <- obj[[8]]$nfam
	nclust <- obj[[8]]$nclust
	dims.use <- obj[[8]]$dims.use
	pos_set <- obj[[8]]$pos_set
	tfl <- obj[[8]]$tfl
	pmat <- obj[[8]]$pmat

    ## set lowlevel fams all to 8
    famind.idr1[famind.idr1>=nfam]=nfam
	
	## reordered tree
    df1 = obj[[10]]
    if(reordertree){
		df1 <- hkm(idrs_all,nclust=nclust,dims.use=dims.use)
#	df1 <- hkm(idrs_all)
    }

	## get the repeats from IDRs we collected
    repeats <- tf.list$fullmatrix[pmat$indices,c(1,33:34)]
    valid_idrs <- which(as.numeric(repeats[,3])>=10) ## the idrs which have a valid repeat => len >=10
	valid_idrs == obj[[9]][[3]]
	

	#posids2 <- posids[which(posids %in%valid)]
	#revposidtable2<-pos_set[which(pos_set[,1] %in%pmat[[2]][posids2]),]
	#revposidtable2<-pos_set[which(pos_set[,1] %in%pmat[[2]][posids_valid]),]  ## this is our gold standard table 
	revposidtable2 <- obj[[3]]

	## positve set ids and idr length > 10  !!! 
	## I need to check if the Repeat is actually the same between the published table and the IDR I have so far

	posids <- obj[[9]][[1]]
	posids_valid <- obj[[9]][[2]]
	cin=c()
	cout=c()
	for(i in 1:dim(repeats[posids_valid,])[1]){
		ind<-which(revposidtable2[,1] %in% repeats[posids_valid,][i,1])
		if(repeats[posids_valid,][i,2] == revposidtable2[ind,6]){
			cin=c(cin,ind)
		}else{
			cout=c(cout,ind)
		}
	}
	revposidtable <-obj[[2]]
#	revposidtable_s <- revposidtable
	reod=c()
	for(i in pmat[[2]][posids]){
		reod=c(reod,which(revposidtable[,1] == i))
	}
#	revposidtable_s <- cbind(revposidtable[reod,],posids,repeats[posids,2:3])
#	revposidtable_s <- cbind(revposidtable_s,revposidtable_s[,6] == revposidtable_s[,8],as.numeric(revposidtable_s[,9]) >=10)
#	colnames(revposidtable_s)[10:11]=c("sameAA","validIDRAA")
#	rownames(revposidtable_s)=1:66
	revposidtable_s <- obj[[6]]

	## finally it is correct
	replist=list();j=0                                                                ## j=  1 2 3   
	for(nn in names(table(revposidtable_s[,6]))){ ## for the repeats in the poslist, nn iselem [A,C,S,...]
		j=j+1;
                replist[[j]]=revposidtable_s[which(revposidtable_s[,6]==nn),7]
	}
	names(replist)=names(table(revposidtable_s[,6]))

	## only the valid ones and the AA from the annotated TF
	replist2=list();j=0                                                                ## j=  1 2 3   
	for(nn in names(table(revposidtable_s[,6]))){ ## for the repeats in the poslist, nn iselem [A,C,S,...]
	    j=j+1;
            replist2[[j]]=revposidtable_s[ which(revposidtable_s[,6]==nn & revposidtable_s[,11] == T),7]
	}
	names(replist2)=names(table(revposidtable_s[,6]))

	## only the valid ones and the AA from the IDR we caught
	onlyv <- which(revposidtable_s[,11] == T)
	replist3=list();j=0                                                                ## j=  1 2 3   
	for(nn in names(table(revposidtable_s[onlyv,8]))){ ## for the repeats in the poslist, nn iselem [A,C,S,...]
	    j=j+1;
            replist3[[j]]=revposidtable_s[ which(revposidtable_s[,8]==nn & revposidtable_s[,11] == T),7]
	}
	names(replist3)=names(table(revposidtable_s[onlyv,8]))


	## now lets make the circle plot colors for the IDRS etc
	mytcol=rep("NA",nsize)
	mytcol2=rep("NA",nsize)
	mytcol3=rep("NA",nsize)
	mytcol4=rep("NA",nsize)  ## colors for repeats we called in our IDRs

	## all I have in my data
	#AAcol=cbind(1:length(names(replist)),names(replist))
	#AAcolF ## has a color foreach of the 20 AA

	allAAin=c()
	namesreplist=names(replist)
	for(j in 1:length(replist)){
	    cj=which(AAcolF[,1] %in%namesreplist[j])
	    mytcol[replist[[j]]]= AAcolF[cj,2]
	}
	allAAin=c(allAAin,namesreplist)
	namesreplist=names(replist2)
	for(j in 1:length(replist2)){
	    cj=which(AAcolF[,1] %in%namesreplist[j])
	    mytcol2[replist2[[j]]]= AAcolF[cj,2]
	}
	allAAin=c(allAAin,namesreplist)
	namesreplist=names(replist3)
	for(j in 1:length(replist3)){
	    cj=which(AAcolF[,1] %in%namesreplist[j])
	    mytcol3[replist3[[j]]]= AAcolF[cj,2]
	}
	allAAin=c(allAAin,namesreplist)
	allAAin=sort(unique(allAAin))

	## get colors for our repeats called in IDRs
	for(j in valid_idrs){
		cj<-which(AAcolF[,1] %in% repeats[j,2])
		if(AAcolF[cj,1] %in% c("A","Q")){
			mytcol4[j]=AAcolF[cj,2]
		}
	}

	## now make are circos plot again
	# df1_labels has the numeric ids in as the leafs in the dendrogram, which are the original indices in pmat

	## this is all for plotting now

	## init plot region
	myth=0.07
	circos.clear()
	#circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.04),gap.degree=0)
	circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0),gap.degree=0)
#	circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0.04),gap.degree=0,canvas.xlim=c(-3,2))
	circos.initialize("foo", xlim = c(0, nsize))

	## depending on the color we supply we get a corresponding circle from outside to inside
	## the colors need to have the same order as the value in the dendrogram is!
	## colors from hkm function using mycolors[1:7] for the the 7 clusters

	## plot TF family has simply base R colors 1 to 8 for the first 8 families,or the color code we give!
	ccol=famind.idr1[df1_labels]
	if(tffam){
	#mtrack(mycol=colsFam[ccol],th=myth)
	mtrack(mycol=colsFam[ccol],th=0.02) ## changed here outer size of tf dbd family
	legend("topleft",paste(    c(tfl$set[1:(nfam-1)],   paste("Rest",sep="")  ), addbraces(table(famind.idr1))   ),pch=15,col=colsFam[1:nfam],cex=cexl,bty="n")
	
	xs=0.75
	text(xs,0.95,"Intrinsically disordered\nregions(IDRs)",adj=0)
	text(xs,1,"Inner circle",adj=0)

	text(xs,0.85,"Homopolymoric amino-\nacid repeats (>10)",adj=0)
	text(xs,0.9,"Middle circle",adj=0)


	text(xs,0.75,"DNA binding domains\n(DBDs)",adj=0)
	text(xs,0.8,"Outer circle",adj=0)
	text(0,0,paste("IDRs\n(",length(df1_labels),")",sep=""))	
	}

	## plot positive set -> colors from mycolors[1:j] where j is the number of different AA repeats found in the data => so A maybe 1, G maybe 2 ... 
	if(1 %in% ppset){
		ccol=mytcol[df1_labels]
		mtrack(mycol=ccol)
	}

	## plot positive set IDRs with paper TF AA
	if(2 %in% ppset){
		ccol=mytcol2[df1_labels]
		mtrack(mycol=ccol)
	}

	## plot positive set IDRs with IDR AA found -> 49 are the same, 4 are different
	if(3 %in% ppset){
		ccol=mytcol3[df1_labels]
		mtrack(mycol=ccol)
	}

	## plot all IRDs with an AA repeat found >= 10
	if(4 %in% ppset){
		ccol=mytcol4[df1_labels]
		mtrack(mycol=ccol,th=myth)
	}

	## this are the clusters
#	ccol=as.numeric(obj$idrs_dall$km$cluster)

	ccol=labels_col(df1)
	mtrack(mycol=colsClust[ccol],th=myth)
	## legend for clusters -> the colors come from hkm function !!!!
	legend("bottomleft",paste(paste("Cluster",1:nclust),  addbraces(table(gene_ids_cluster))),pch=16,col=colors7[1:nclust],cex=cexl,bty="n")
	#legend("bottomleft",paste(paste("Cluster",1:nclust),  addbraces(table(gene_ids_cluster))),pch=16,col=colsClust[1:nclust],cex=cexl,bty="n")

	## legend for tf famlies

	to=c(); ## should hold then increasing indices anyways
	for(i in allAAin){ 
		to=c(to,which(AAcolF[,1] %in% i))    
	}

	## make the legend matrix
	nmat=matrix(nrow=length(allAAin),ncol=3)
	rownames(nmat)=allAAin
	for(i in 1:length(allAAin)){
		nn=allAAin[i]
		ind <- which(names(replist)==nn)
		if(length(ind) > 0){
			nmat[i,1]=length(replist[[ind]])
		}else{
			nmat[i,1]=0
		}
	}

	for(i in 1:length(allAAin)){
		nn=allAAin[i]
		ind <- which(names(replist2)==nn)
		if(length(ind) > 0){
			nmat[i,2]=length(replist2[[ind]])
		}else{
			nmat[i,2]=0
		}
	}

	for(i in 1:length(allAAin)){
		nn=allAAin[i]
		ind <- which(names(replist3)==nn)
		if(length(ind) > 0){
			nmat[i,3]=length(replist3[[ind]])
		}else{
			nmat[i,3]=0
		}
	}

	fnl=c()
	for(i in 1:dim(nmat)[1]){
		fnl=c(fnl,paste(rownames(nmat)[i],paste0(nmat[i,],collapse=" ")))
	}

	#legend("bottomleft",fnl,col=AAcolF[to,2],lwd=2,cex=cexl)
	#table(repeats[valid_idrs,2])
	rti<-table(repeats[valid_idrs,2])
	fnl=paste("Poly",names(rti)," ",addbraces(as.numeric(rti)),sep="")
	fnlc=c()
	mykeep=c()
	jo=0
	for(j in names(rti)){
		cj=which(AAcolF[,1] %in% j)
		fnlc=c(fnlc,AAcolF[cj,2])
		jo=jo+1
		if(AAcolF[cj,1]%in%c("A","Q")){
			mykeep=c(mykeep,jo)
		}
	}


	## add foxp
	dend=hkm(obj[[7]],nclust=nclust,dims.use=dims.use)

	if(hlf==T){
		mytcol=rep("white",nsize)
		FOXP=c();for(ii in 1080:1083){
		FOXP=c(FOXP,which(labels(dend) %in% ii))}
		mytcol[FOXP[3]]="red"
		mytcol[FOXP[c(1,2,4)]]="blue"
#		mtrack(mycol=mytcol,th=0.04)
		mytcol[FOXP[c(1,2)]]="cornflowerblue"
		mtrack(mycol=mytcol,th=0.04)
	}
	if(length(hlfg)>0){
		mytcol=rep("white",nsize)
		gn<-names(obj$idrs_dall$km$cluster) ## original order
		gni = which(gn %in% hlfg) ## get index in this
		posi=which(df1_labels %in% gni) ## get now the position in dendrogram
		mytcol[posi]="darkolivegreen"
		mtrack(mycol=mytcol,th=0.04)
	}


	if(4 %in% ppset){		
		legend("bottomright",fnl[mykeep],col=fnlc[mykeep],lwd=2,cex=cexl,bty="n")
	}
	if(plotdend){
		 max_height = attr(dend, "height")
		circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
		circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
        circos.dendrogram(dend, max_height = max_height)}, track.height = 0.5, bg.border = NA)
		#circos.dendrogram(dend, labels_track_height = NA, dend_track_height = .3)
	}
	return(gene_ids_cluster)
}


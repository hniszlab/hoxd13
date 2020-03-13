## this is the length of a list of genes 
#gene.lengths=length(c(gl,gl2))

## number of clusters 
#clustl=11

## make dot plot as in bla
getmat <-function(p=p.obj,g=c(gl,gl2),clustl=11,mfun="sum"){
	gene.lengths=length(g)
	## allocate matrix
	dfm <- matrix(nrow=clustl*gene.lengths,ncol=6)
	## set names
	colnames(dfm)=c("Gene","Cluster","Exprs","Excell","Frac","Totalc")
	j=0
	for(gene in g){
		j=j+1
		gene.i <- which(rownames(p@data) == gene)
		if(length(gene.i) > 0){
			cpef <-list()
			for(i in 1:11){
				cpef[[i]] <- as.vector(p@data[gene.i,p@ident==i])
			}
			names(cpef)=paste("C",1:11,sep="")
			x=cpef
			nlabels <- unlist(lapply(x,length))
			isex <- nlabels-unlist(lapply(lapply(x,is0),length))

			## which of those are expressed
			islex<-list()
			medex=c()

			## what function to apply ## name from mfun and needs to match R function
			ffun <- get_var(mfun)

			for(i in 1:11){
				islex[i]=list(which(x[[i]] > 0))
				#get_median_per cluster of expressed cells
				#medex=c(medex,median(x[[i]][islex[[i]]]))

				## get sum of expressed reads now 
				medex=c(medex,ffun(x[[i]][islex[[i]]]))
			}
			## normalize now by highest expression values
			medex=medex/max(medex)


			from=((j-1)*clustl)+1
			to=j*clustl
			dfm[(from:to),1]=j
			dfm[(from:to),2]=1:clustl
			dfm[(from:to),3]=medex
			dfm[(from:to),4]=unlist(lapply(islex,length))
			dfm[(from:to),6]=unlist(lapply(x,length))
			dfm[(from:to),5]=dfm[(from:to),4]/dfm[(from:to),6]
		}else{
			print(paste(gene,"not in matrix"))
			from=((j-1)*clustl)+1
			to=j*clustl
			dfm[(from:to),1]=j
			dfm[(from:to),2]=1:clustl
			dfm[(from:to),3]=0
			dfm[(from:to),4]=0
			dfm[(from:to),6]=unlist(lapply(x,length))
			dfm[(from:to),5]=0

		}

		#dfm[(from:to),5]=rep(gene,clustl)
	}
	dfm[is.na(dfm)]=0
	dfmdf <- data.frame(dfm)
	return(dfmdf)
}


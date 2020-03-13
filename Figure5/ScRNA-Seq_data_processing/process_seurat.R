process <- function(obj){
		## rows is a 3 column matrix with gene symbol, ensemblID and refseq
	    rows.f <- rownames(obj@raw.data)
		## get ens_gene_ids that are in the mitochondrial gene fraction
	    rows.mito.i <- which(rows[,1]%in%mito)
	    ## table with ensgene and gene symbol only
	    rows.mito<-rows[rows.mito.i,1:2]
	    ## indices in rows all that are still in after filtering
	    mito.genes <- which(rows.f%in%rows.mito[,2])
		percent.mito <- Matrix::colSums(obj@raw.data[mito.genes, ])/Matrix::colSums(obj@raw.data)

		obj <- AddMetaData(object = obj, metadata = percent.mito, col.name = "percent.mito")
		## take out cells with less than 200 genes or more then 7000 and more than 5% of mitochondrial reads
		obj <- FilterCells(object = obj, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(7000, 0.05))
		## normalize now
		obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
		
		## variable genes
		obj <- FindVariableGenes(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

		## scale
		obj            <- ScaleData(object = obj)
		## pca
		obj            <- RunPCA(object = obj, pc.genes = obj@var.genes)
		## tsne
		obj            <- RunTSNE(object = obj, dims.use = 1:10, do.fast = TRUE,reduction.use = "pca")
		## clusters
		obj            <- FindClusters(object = obj, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
		return(obj)
}

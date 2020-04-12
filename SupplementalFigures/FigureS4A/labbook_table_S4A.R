#!/usr/bin/Rscript

## start an R session and run 
#../../Figure5/ScRNA-Seq_data_processing/generate_figure_5B_and_analyses.R
## before continuing with the code below


## get all genes expressed in cluster 1:11 with mean expression > 0,0.5 and 1
for(thres in c(0,0.5,1)){
	allg <- c()
	for(i in 0:10){
		clustc <- which(pbmc.ctrl@ident == i)
		allg <-c(allg,length(which(rowMeans(pbmc.ctrl@data[,clustc]) > thres)))
	}
	print(allg)
}

table(pbmc.ctrl@ident)
# 0   1   2   3   4   5   6   7   8   9  10
#722 614 593 564 564 540 410 301  70  59  27

round(100*table(pbmc.ctrl@ident)/4464,1)

#   0    1    2    3    4    5    6    7    8    9   10 
# 16.2 13.8 13.3 12.6 12.6 12.1  9.2  6.7  1.6  1.3  0.6 

#>0: 15719 15555 15931 15372 15415 15224 15045 13884 13239 12695  8618
#>0.5 2046 2212 2163 1913 1900 1781 1881 1545 2013 1644  452
#>1     823 864 810 717 745 650 731 631 805 580 194



#making Plot inset --------------------------------------------------------------

install.packages(directlabels)
require(ggplot2)
require(directlabels)


#making Plot p-------------------------------------------------------------------

#setwd
setwd("/Users/basu/Desktop/TF\ Phase\ Final\ Submission/Figure\ 3/Figure\ 3A-D/3B-C_Reanalysis_200306")
#specify path to folder containing analysis files


#readTable
red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("D13WT-2uM", red[i,5])==TRUE)
  {red[i,6]="Hoxd13WT 2 uM"}
  else if(grepl("D137A", red[i,5])==TRUE)
  {red[i,6]="Hoxd137A 2 uM"}
  else if(grepl("D1310A", red[i,5])==TRUE)
  {red[i,6]="Hoxd1310A 2 uM"}
  else{red[i,6]="Hoxd13WT 4 uM"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix, levels = c("Hoxd1310A 2 uM", "Hoxd137A 2 uM","Hoxd13WT 2 uM","Hoxd13WT 4 uM"))

Hoxd13 <- red

Hoxd13$mix <- as.factor(Hoxd13$mix)

Hoxd13$mix <- factor(Hoxd13$mix, levels = c("Hoxd1310A 2 uM", "Hoxd137A 2 uM","Hoxd13WT 2 uM","Hoxd13WT 4 uM"))


#make plot 3B
p <- ggplot(Hoxd13, aes(y=red, x=green, color = Hoxd13$mix)) + geom_point(alpha=0.3, size = Hoxd13$area/500000) + theme_classic() + scale_color_manual(values=c('red', 'orange', 'lightseagreen', 'blue'))

#make plot 3C
p2 <- ggplot(Hoxd13, aes(x=mix, y=green/red)) + theme_classic() + geom_jitter() + ylab('Med1 content / D13 content')




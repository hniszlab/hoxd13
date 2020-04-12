#making Plot inset --------------------------------------------------------------

install.packages(directlabels)
require(ggplot2)
require(directlabels)


#making Plot p-------------------------------------------------------------------

#set wd
setwd()
#specify path to folder containing analysis files

#read table
red <- read.csv2("objectsGREEN", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("D13WT-2uM", red[i,5])==TRUE)
  {red[i,6]="Hoxd13WT 2 uM"}
  else if(grepl("D137A", red[i,5])==TRUE)
  {red[i,6]="Hoxd137A 2 uM"}
  else if(grepl("D1310A", red[i,5])==TRUE)
  {red[i,6]="Hoxd1310A 2 uM"}
  else if(grepl("mCherry", red[i,5])==TRUE)
  {red[i,6]="mCherry 2 uM"}
  else{red[i,6]="Hoxd13WT 4 uM"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix, levels = c("mCherry 2 uM","Hoxd13WT 2 uM","Hoxd13WT 4 uM", "Hoxd137A 2 uM","Hoxd1310A 2 uM"))

Hoxd13 <- red

Hoxd13$mix <- as.factor(Hoxd13$mix)

Hoxd13$mix <- factor(Hoxd13$mix, levels = c("mCherry 2 uM","Hoxd13WT 2 uM","Hoxd13WT 4 uM", "Hoxd137A 2 uM","Hoxd1310A 2 uM"))

p <- ggplot(Hoxd13, aes(y=red, x=green, color = Hoxd13$mix)) + geom_point(alpha=0.3, size = Hoxd13$area/500000) + theme_classic() + scale_color_manual(values=c('red', 'orange', 'lightseagreen'   , 'blue'))

#make Fig 3D
p2 <- ggplot(Hoxd13, aes(x=mix, y=red)) + geom_boxplot(outlier.shape = NA) + ylab('mCherry') + ylim(-2,260) + theme_classic()

p3 <- ggplot(Hoxd13, aes(x=mix, y=red)) + geom_jitter() + ylab('mCherry') + theme_classic()


# normalized DFP
  
normHoxd13 <- Hoxd13

normHoxd13$normG <- normHoxd13$green/max(Hoxd13$green)
normHoxd13$normR <- normHoxd13$red/max(Hoxd13$red)

n <- ggplot(normHoxd13, aes(y=normR, x=normG, color = Hoxd13$mix)) + geom_point(alpha=0.3, size = Hoxd13$area/500000) + theme_classic() + scale_color_manual(values=c('red', 'orange', 'lightseagreen', 'blue'))



View(Hoxd13)

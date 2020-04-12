#setwd to folder containing 6H Image analysis files
setwd("")

#read table of values
red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("A13_3,5uM", red[i,5])==TRUE)
  {red[i,6]="Hoxa13WT 3,5 uM"}
  else{red[i,6]="Hoxa13(+7A) 3,5uM"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix)

#make Figure 6H Dual Flourescence Plot
ggplot(red, aes(y=red, x=green, color=mix)) + geom_point(alpha=0.4, size = red$area/500000) + theme_classic() + ylab('HOXD13 condensates [mCherry signal]') + xlab('Mediator enrichment [GFP signal]') + scale_color_manual(values=c('red','blue'))






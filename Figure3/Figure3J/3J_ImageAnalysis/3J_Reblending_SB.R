#setwd 
setwd()
#specify path to folder containing analysis files

#read data
red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("Mito100", red[i,5])==TRUE)
  {red[i,6]="100 uM Mito HCl"}
  else if(grepl("LipoAcid", red[i,5])==TRUE)
  {red[i,6]="100 uM Lipoic Acid"}
  else if(grepl("LipoAmide", red[i,5])==TRUE)
  {red[i,6]="100 uM Lipoamide"}
  else if(grepl("D13WT_2uM_Med1_1uM_Vehicle", red[i,5])==TRUE)
  {red[i,6]="WT DMSO"}
  else if(grepl("D137A_2uM_Med1_1uM_Vehicle", red[i,5])==TRUE)
  {red[i,6]="7A DMSO"}
  else{red[i,6]="DMSO"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix, levels = c("100 uM Mito HCl", "100 uM Lipoic Acid", "100 uM Lipoamide", "7A DMSO", "WT DMSO"))
ggplot(red, aes(y=red, x=green, color=mix)) + geom_point(alpha=0.4, size = red$area/500000) + theme_classic() + ylab('HOXD13 condensates [mCherry signal]') + xlab('Mediator enrichment [GFP signal]') + scale_color_manual(values=c('red','lightseagreen','purple', 'orange','blue'))

red2 <- subset(red, green/red <  8)

#make Fig 3J
ggplot(red2, aes(x=mix, y=green/red)) + theme_classic() + geom_jitter() + ylab('Med1 content / D13 content') + scale_y_log10()
 

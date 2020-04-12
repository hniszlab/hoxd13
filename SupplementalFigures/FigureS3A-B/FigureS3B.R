#main -------------------------------
require(ggplot2)
require(ggthemes)

#setwd to FigS3B data folder
setwd()

#read data
red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("wt_AF568", red[i,5])==TRUE)
  {red[i,6]="WT"}
  else if(grepl("wt_HD6%", red[i,5])==TRUE)
  {red[i,6]="WT HD6%"}
  else if(grepl("homozygous_AF568", red[i,5])==TRUE)
  {red[i,6]="SPDH"}
  else if(grepl("homozygous_HD6%", red[i,5])==TRUE)
  {red[i,6]="SPDH HD6%"}
  else if(grepl("spdh_", red[i,5])==TRUE)
  {red[i,6]="SPDH"}
  else if(grepl("spdh_hexanediol_", red[i,5])==TRUE)
  {red[i,6]="SPDH HD6%"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix, levels = c("WT", "WT HD6%", "SPDH", "SPDH HD6%"))
redmain <- subset(red, red$mix == "WT" | red$mix == "SPDH")


#Make Fig S3B
p <- ggplot(redmain, aes(y = redmain$Channel.1, x = redmain$mix)) + geom_boxplot(outlier.shape = NA) + geom_rangeframe() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + 
  xlab('') + ylab('Hoxd13 Puncta Intensity [AU]') 


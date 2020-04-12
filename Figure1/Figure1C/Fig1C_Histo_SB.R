#set wd, to folder "Fig1C_VutaraSRX_ClusterAnalysis" 
setwd('')

#read data
gatedcell <- read.csv2("Fig1C-TableofValues.csv", header=T, sep = ",", dec = ".")

cell2 <- subset(gatedcell, gatedcell$ID == 'cell2')

#Fig 1C Inset
require(ggplot2)
require(RColorBrewer)

ggplot(cell2, aes(x=size, color = ID, fill = ID)) + theme_classic() + xlab("2 x Radius of Gyration (nm)") +
  geom_density(alpha=.4) + scale_fill_manual(values=c('purple','blue')) + scale_color_manual(values=c('purple','blue')) +
  ggtitle("STORM : Hoxd13 Condensate Size  \nParticle Count 15 - 50")

summary(cell2$size)



#setwd, specify path to folder containing analysis files ('4F_ImageAnalysis')

setwd()


#read values

merge3 <- read.csv("Hoxd13tether_LacIO_TotalMerge_Gated_200115.csv", sep = ",")
merge4 <- subset(merge3, merge3$mix == "Med1" | merge3$mix == "Hoxd13 WT" | merge3$mix == "Hoxd13 7A" | merge3$mix == "Hoxd13 10A")


#Make Fig 4F

ggplot(merge4, aes(x=mix, y=(merge4$TetherYellow/merge4$RingYellow))) + theme_classic() + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) + ylab('yellow enrichment at tether') + ylim(0.9,1.6) 


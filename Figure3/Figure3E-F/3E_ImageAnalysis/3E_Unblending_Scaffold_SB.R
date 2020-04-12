#setwd
setwd()
#specify path to folder containing analysis files

install.packages(directlabels)
require(ggplot2)
require(directlabels)

#make insets --------------------------------------------------------------

z <- read.csv2('Hoxd13_PhaseDiagram_summary', dec = '.', head = T, sep = ",")

g <- subset(z, z$Mean.Intensity > 40 & z$Mean.Intensity < 80)

t <- subset(z, z$ID == '1 Hoxd13(WT)' | z$ID == '1 Hoxd13(10A)')
L <- rbind(g,t)

inset1 <- ggplot(z, aes(z$Concentration..uM., z$Mean.Intensity,color=z$Genotype)) + geom_line() + theme_classic() +
  geom_hline(yintercept = mean(g$Mean.Intensity), linetype='dotted') + ylab('AU') + xlab('[uM]') +
  geom_vline(xintercept = 1, linetype='dotted') +
  labs(color = 'Genotype') + 
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'Red', 'Orange', 'Blue')) +
  geom_point(data=L, aes(x=L$Concentration..uM., y=L$Mean.Intensity, color = L$ID), size = 4, shape = 0) + 
  geom_point(data=L, aes(x=L$Concentration..uM., y=L$Mean.Intensity, color = L$ID)) + 
  theme(legend.position='none') + ggtitle('Hoxd13 Phase Diagram') + 
  geom_dl(aes(label = z$Genotype), method = list(cex = 0.5, dl.trans(y = y + 0.2, x = x - 1.1),"last.points")) + ylim(0,220)


#making Plot p-------------------------------------------------------------------

red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("Hoxd13-wt-n-1uM_Med1", red[i,5])==TRUE)
  {red[i,6]="Hoxd13WT 1 uM"}
  else if(grepl("Hoxd13-wt-n-5uM_Med1-1uM_10PEG_preas2_mix2", red[i,5])==TRUE)
  {red[i,6]="Hoxd13WT 5 uM"}
  else if(grepl("Hoxd137A-1uM_Med1-1uM", red[i,5])==TRUE)
  {red[i,6]="Hoxd137A 1 uM"}
  else if(grepl("Hoxd1310A-1uM_Med1-1uM", red[i,5])==TRUE)
  {red[i,6]="Hoxd1310A 1 uM"}
  else if(grepl("Hoxd1310A-0-2uM_Med1", red[i,5])==TRUE)
  {red[i,6]="Hoxd1310A 0.2 uM"}
  else{red[i,6]="Misc"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix, levels = c("Hoxd13WT 1 uM","Hoxd13WT 5 uM","Hoxd137A 1 uM","Hoxd1310A 0.2 uM","Hoxd1310A 1 uM",   "Misc"))

Hoxd13 <- subset(red,mix == "Hoxd13WT 1 uM"| mix == "Hoxd13WT 5 uM" | mix == "Hoxd137A 1 uM" | mix == "Hoxd1310A 0.2 uM" | mix == "Hoxd1310A 1 uM")

p <- ggplot(Hoxd13, aes(y=red, x=green)) + geom_point(alpha=Hoxd13$green/256, size = Hoxd13$area/500000, color = 'green') + 
  geom_point(size = Hoxd13$area/500000, fill = NA, shape = 1) +
  theme_classic() + ylab('HOXD13 condensates [mCherry signal]') + xlab('Mediator enrichment [GFP signal]') + 
  scale_x_log10() + annotation_logticks(sides = 'b')  + ylim(0,150) + labs(color = '1 uM Med1 +') + facet_wrap(~mix,ncol=1) +
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd1310A 0.2 uM',], aes(xintercept=median(green)), linetype = 'dashed') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd137A 1 uM',], aes(xintercept=median(green)), linetype = 'dashed') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd13WT 5 uM',], aes(xintercept=median(green)), linetype = 'dashed') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd13WT 1 uM',], aes(xintercept=median(green)), linetype = 'dashed') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd1310A 1 uM',], aes(xintercept=median(green)), linetype = 'dashed') +
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd1310A 0.2 uM',], aes(xintercept=mean(green)), linetype = 'dotted') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd137A 1 uM',], aes(xintercept=mean(green)), linetype = 'dotted') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd13WT 5 uM',], aes(xintercept=mean(green)), linetype = 'dotted') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd13WT 1 uM',], aes(xintercept=mean(green)), linetype = 'dotted') + 
  geom_vline(data=Hoxd13[Hoxd13$mix=='Hoxd1310A 1 uM',], aes(xintercept=mean(green)), linetype = 'dotted') 


#make Fig3E---------------------------

R = p + annotation_custom(grob=ggplotGrob(inset1 + 
  theme(plot.title = element_text(size = 6),axis.text=element_text(size=6),axis.title=element_text(size=6), axis.text.y=element_blank(), panel.background = element_rect(fill = "transparent", colour = NA),plot.background = element_rect(fill = "transparent", colour = NA))), ymin = 25, ymax = 160, xmin = -1, xmax = 0.1) + 
  ggtitle("Preassembled Med1 [1 uM] +")





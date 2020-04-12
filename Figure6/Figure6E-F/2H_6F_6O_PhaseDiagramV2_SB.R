#Compile the three functions, ReadTable, PhaseDiagramDiameter, PhaseDiagramIntensity

#ReadTable <- function(gene)
MasterTable <- function(gene)
{
  #setwd path to '/PhaseDiagram_shaonR', or folder containing Runx2, Hoxa13, Hoxd13 subdirectories
  setwd(paste('',gene,sep=""))
  require(ggplot2)
  require(RColorBrewer)
  condensates2 <- read.csv2("PhaseDiagram_Values", sep = ",", dec = ".", header=T)
  return(condensates2)
}

PhaseDiagramDiameter <- function(gene,y)
{
  if(gene =='Hoxa13')
  {
    y <- y[!(y$ID=='0,01 Hoxa13(7A)' & y$Diameter > 1000),]
  }
  
  z <- aggregate(y[,c("diameteradjusted","Concentration..uM.")],list(y$ID), mean)
  
  sd <- aggregate(y[,c("diameteradjusted","Concentration..uM.")],list(y$ID), sd)
  
  int <- aggregate(y[,c("Intensity", "Concentration..uM.")],list(y$ID), mean)
  
  for(i in 1: nrow(z))
  {
    value <- z[i,1]
    if(grepl("WT", value)==TRUE)
    {z[i,4]=paste(gene,"(WT)",sep="")}
    else if(grepl("7A", value)==TRUE)
    {z[i,4]=paste(gene,"(7A)",sep="")}
    else if(grepl("10A", value)==TRUE)
    {z[i,4]=paste(gene,"(10A)",sep="")}
    else
    {z[i,4]=paste(gene,"(WT)",sep="")}
  }
  
  z$errorminus <- z$diameteradjusted - sd$diameteradjusted
  
  z$errorplus <- z$diameteradjusted + sd$diameteradjusted
  
  if(gene == 'Hoxd13')
  {
    z <- subset(z, z$Concentration..uM.!=10)
    int <- subset(int, int$Concentration..uM.!=10)
    color <- c('red', 'orange', 'blue')
    touch = 0.7
  }
  else if(gene == 'Runx2')
  {
    z <- subset(z, z$Concentration..uM.!=1)
    int <- subset(int, int$Concentration..uM.!=1)
    color <- c('red', 'blue')
    touch = 0.7
  }
  else if(gene == 'Hoxa13')
  {
    color <- c('red', 'blue')
    touch = 0.8
  }
  
  write.csv(z, paste(gene,"MicroscopyReduxDiameter",sep=""))
  
  p <- ggplot(z, aes(x=as.factor(z$Concentration..uM.), y = z$V4, color = as.factor(z$V4))) + 
    geom_point(alpha = touch*(int$Intensity/(max(int$Intensity))), size = (z$errorplus/100)) + 
    geom_point(alpha = touch*(int$Intensity/(max(int$Intensity))), size = (z$diameteradjusted/100)) + 
    geom_point(alpha = touch*(int$Intensity/(max(int$Intensity))), size = (z$errorminus/100)) +
    ylab(NULL) + xlab(NULL) + theme_classic() + scale_color_manual(values = color) +
    #ggtitle (paste(gene,': Phase Diagram (Condensate Diameter)', sep=""))
    ggtitle('Mean Projection')
  
  #save pdf 
  ggsave(p, filename = paste(gene, "_Diameter.pdf", sep = ""))
  return(p)
}

PhaseDiagramIntensity = function(gene, y)
{

  if(gene == 'Hoxd13')
  {
    y <- subset(y, y$Concentration..uM.!=10)
    color <- c('red', 'orange', 'blue')
  }
  else if(gene == 'Runx2')
  {
    y <- subset(y, y$Concentration..uM.!=1)
    color <- c('red', 'blue')
  }
  else if(gene == 'Hoxa13')
  {
    color <- c('red', 'blue')
  }
  
p <- ggplot(y, aes(as.factor(y$Concentration..uM.), y$Intensity, color = y$genotype))  +
  geom_point(position = position_jitterdodge(), alpha = 0.4, size = y$diameteradjusted/1000) + theme_classic() +
  geom_boxplot(outlier.size =0, fill = NA) + scale_color_manual(values=color) + ylab('Flour Intensity [mCherry]') +
  xlab('Concentration [uM]') + ggtitle(paste(gene,': Phase Diagram', sep="")) + labs(color = 'Genotype') + theme(plot.title = element_text(size=18))
}



#Misc ------------------------------------------------------- after all else compiled

Boss <- function(gene)
#options  are 'Hoxd13', 'Hoxa13', or 'Runx2'
#this function will return output for Fig 3H, 6F, or 6O
{
  y <- MasterTable(gene)
  z <- PhaseDiagramDiameter(gene,y)
  p <- PhaseDiagramIntensity(gene,y)
  if(gene == 'Hoxd13')
  {
  L = p + annotation_custom(grob=ggplotGrob(z + theme(legend.position='none', axis.text.y=element_blank(), plot.title = element_text(size = 10),
                                                      axis.ticks.y = element_blank())), 
                            ymax = 265, ymin = 180, xmin = 0.5, xmax = 1.6)
  }
  else
  {
  L = p + annotation_custom(grob=ggplotGrob(z + theme(legend.position='none', axis.text.y=element_blank(), plot.title = element_text(size = 10),
                                                        axis.ticks.y = element_blank())), 
                              ymin = max(y$Intensity)-(max(y$Intensity)/3), ymax = max(y$Intensity)+10, xmin = 0.5, xmax = 3)
  }
  return(L)
}


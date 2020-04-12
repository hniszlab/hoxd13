
#setwd to folder Fig1J_ImageAnalysis
setwd('')



#read data------------------------------------------------


z <- read.csv2("summary", header=T, sep = ";", dec=",")

for(i in 1: nrow(z))
{
  value <- z[i,c("name")]
  if(grepl("11", value)==TRUE)
  {z[i,c("ID")]='Hoxd13 12 PEG 10'}
  else
  {z[i,c("ID")]='Hoxd13 6 PEG 10'}
}

phase <- aggregate(z[,c("Area..Area..R")], list(z$ID), mean)
colnames(phase) <- c("ID", "area")

for(i in 1: nrow(phase))
{
  value <- phase[i,1]
  if(grepl("12", value)==TRUE)
  {phase[i,c("concentration (uM)")]=12
  phase[i,c("Crowding Agent (%PEG)")]=10
  phase[i,c("condensates")]='yes'}
  else
  {phase[i,c("concentration (uM)")]=6
  phase[i,c("Crowding Agent (%PEG)")]=10
  phase[i,c("condensates")]='yes'}
}

neg <- data.frame("Hoxd13 12 PEG 5",0,12,5,"no")
colnames(neg) <- colnames(phase)
neg2 <- data.frame("Hoxd13 6 PEG 5",0,6,5,"no")
colnames(neg2) <- colnames(phase)
neg3 <- data.frame("Hoxd13 2 PEG 5",0,2,5,"no")
colnames(neg3) <- colnames(phase)
pos <- data.frame("Hoxd13 2 PEG 10",0,2,10,"yes")
colnames(pos) <- colnames(phase)

newdf1 <- rbind(phase,pos,neg,neg2,neg3)

newdf1$`concentration (uM)` <- as.factor(newdf1$`concentration (uM)`)
newdf1$`Crowding Agent (%PEG)` <- as.factor(newdf1$`Crowding Agent (%PEG)`)
newdf2 <- newdf1

#Make Fig 1J------------------------------------------------
ggplot(newdf2, aes(y=`Crowding Agent (%PEG)`, x=`concentration (uM)`,color=condensates)) + geom_point(alpha=0.8, size = 2 + newdf2$area/35000) + theme_classic() + scale_color_manual(values =c('grey',"purple")) + xlab("concentration (uM)") + ylab("Crowding Agent (% PEG)")


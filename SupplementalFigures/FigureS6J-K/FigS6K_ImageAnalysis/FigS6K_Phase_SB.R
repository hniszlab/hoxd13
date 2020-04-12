
#setwd to folder S6K IMAGE ANALYSIS, containing values.csv
setwd()

#read data------------------------------------------------

z <- read.csv2("values", header=T, sep = ",", dec=",")

for(i in 1: nrow(z))
{
  value <- z[i,c("name")]
  if(grepl("40uM_TBP_10%PEG", value)==TRUE)
  {z[i,c("ID")]='TBP 40 PEG 10'}
  else
  {z[i,c("ID")]='TBP 20 PEG 10'}
}

phase <- aggregate(z[,c("Area..Area..R")], list(z$ID), mean)
colnames(phase) <- c("ID", "area")

for(i in 1: nrow(phase))
{
  value <- phase[i,1]
  if(grepl("40", value)==TRUE)
  {phase[i,c("concentration (uM)")]=40
  phase[i,c("Crowding Agent (%PEG)")]=10
  phase[i,c("condensates")]='yes'}
  else
  {phase[i,c("concentration (uM)")]=20
  phase[i,c("Crowding Agent (%PEG)")]=10
  phase[i,c("condensates")]='yes'}
}

neg <- data.frame("TBP 40 PEG 5",0,40,5,"no")
colnames(neg) <- colnames(phase)
neg2 <- data.frame("TBP 20 PEG 5",0,20,5,"no")
colnames(neg2) <- colnames(phase)
neg3 <- data.frame("TBP 10 PEG 5",0,10,5,"no")
colnames(neg3) <- colnames(phase)
neg4 <- data.frame("TBP 10 PEG 10",0,10,10,"no")
colnames(neg4) <- colnames(phase)

newdf1 <- rbind(phase,neg,neg2,neg3,neg4)

newdf1$`concentration (uM)` <- as.factor(newdf1$`concentration (uM)`)
newdf1$`Crowding Agent (%PEG)` <- as.factor(newdf1$`Crowding Agent (%PEG)`)
newdf2 <- newdf1

#Make Fig S2K------------------------------------------------

ggplot(newdf2, aes(y=`Crowding Agent (%PEG)`, x=`concentration (uM)`,color=condensates)) + geom_point(alpha=0.8, size = 2 + newdf2$area/35000) + theme_classic() + scale_color_manual(values =c('grey',"purple")) + xlab("concentration (uM)") + ylab("Crowding Agent (% PEG)")


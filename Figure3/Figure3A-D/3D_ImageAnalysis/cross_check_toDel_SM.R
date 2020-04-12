
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

summary(Hoxd13)

mc=which(Hoxd13$mix == "mCherry 2 uM")
wt1=which(Hoxd13$mix == "Hoxd13WT 2 uM")
wt2=which(Hoxd13$mix == "Hoxd13WT 4 uM")
h10a=which(Hoxd13$mix == "Hoxd1310A 2 uM")
h7a=which(Hoxd13$mix == "Hoxd137A 2 uM")

lc=list(Hoxd13$red[mc],Hoxd13$red[wt1],Hoxd13$red[wt2],Hoxd13$red[h7a],Hoxd13$red[h10a])


pdf("Figure3D_repro_SM_toDelete.pdf",12,6)
par(mfcol=c(1,2))
plot(ecdf(Hoxd13$red[wt1]),verticals=T,pch=".",xlab="mCherry fluorescence intensity (U)",main="MED1-IDR-GFP-containing condensates",ylab="Cumulative PDF",las=1,xlim=c(0,400))
lines(ecdf(Hoxd13$red[h7a]),col="red",lwd=1,verticals=T,pch=".")
lines(ecdf(Hoxd13$red[h10a]),col="brown",lwd=1,verticals=T,pch=".")
lines(ecdf(Hoxd13$red[wt2]),col="grey",lwd=1,verticals=T,pch=".")
lines(ecdf(Hoxd13$red[mc]),col="orange",lwd=1,verticals=T,pch=".")
nn=c("mc","wt1","wt2","7A","10A")
nn2= paste(nn,unlist(lapply(lc,length)),sep=", n=")
legend("bottomright",nn2,col=c("orange","black","grey","red","brown"),lwd=1)

boxplot(lc,col=c("orange","black","grey","red","brown"),names=c("mc","wt1","wt2","7A","10A"))
dev.off()

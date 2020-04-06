#!/usr/bin/Rscript

#in ucsc browser select region to get values for
#-> table browser
#-> merged track
#-> get values

#run thing below
# run get values script

gene="Hoxd12"
ymax=115

all=list()
j=0
for(clust in c("ctrl","spdh")){
	j=j+1
	all[[j]]=list()
	for(i in 1:11){
		all[[j]][[i]]=read.table(paste("cluster_",clust,"_",i,"_BPM_",gene,".wig",sep=""))
	}
	i=12
	all[[j]][[i]]=read.table(paste("cluster_",clust,"_merged_BPM_",gene,".wig",sep=""))
}

## add rows in missing clusters
allcr<-list()

for(j in 1:2){
allcr[[j]]=list()
for(i in 1:11){
	mn=all[[1]][[12]] ## this is our reference

	## we make new list here now

	allcr[[j]][[i]]=mn
	allcr[[j]][[i]][,4]=0 ## orig value from merged, but col4 is 0
	
	
	inc <- which(all[[j]][[i]][,2] %in%mn[,2])
	if(length(inc) >0){
		tmp=all[[j]][[i]][inc,]
		myind <- which(mn[,2]%in%tmp[,2])
	## add orig values from orig matrix now
		allcr[[j]][[i]][myind,4]=tmp[,4]
	}

}
allcr[[j]][[12]]=all[[j]][[12]]
}

#svglite(paste(gene,"_replot.svg",sep=""),6,18)
pdf(paste(gene,"_replot.pdf",sep=""),6,18)

num1=c()
num2=c()
num3=c()
num4=c()

par(mfcol=c(7,2),oma=c(0.1,2,5,2),mar=c(0.1,2,5,2))
for(i in c(12,1:11)){
	v1=mean(all[[1]][[i]][,4])
	m1=max(all[[1]][[i]][,4])

	v2=mean(all[[2]][[i]][,4])
	m2=max(all[[2]][[i]][,4])
	
	num1=c(num1,m1)
	num2=c(num2,m2)
	num3=c(num3,v1)
	num4=c(num4,v2)

	if(i==12){
		namex="merged"
	}else{
		namex=paste("cluster",i)
	}

 barplot(allcr[[1]][[i]][,4],las=1,ylim=c(0,ymax),col=rgb(0,64,255,maxColorValue=256),space=c(0,0),border=rgb(0,64,255,maxColorValue=256),main=paste(namex,"\nmean ctrl=",round(v1),"spdh=",round(v2),"\nmax ctrl=",m1,"spdh=",m2))
 barplot(allcr[[2]][[i]][,4],add=T,col=rgb(255,64,0,maxColorValue=256),ylab="",yaxt="n",space=c(0,0),border=rgb(255,64,0,maxColorValue=256))
	if(v1 > v2){
		barplot(allcr[[2]][[i]][,4],add=T,col="black",ylab="",yaxt="n",space=c(0,0))
		}else{
		barplot(allcr[[1]][[i]][,4],add=T,col="black",ylab="",yaxt="n",space=c(0,0))
		}


}
#mvec=paste(round(num1),round(num2),sep=" / ")
#mvec2=c("merged",paste("cluster",1:11))

mvecwt0=round(num1)
mvecspdh0=round(num2)

mvecWT=round(num3)
mvecSPDH=round(num4)
mvec2=c("merged",paste("cluster",1:11))

plot(1,xlim=c(1,100),ylim=c(1,130),col="white",main="mean values of regions",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
for(i in c(12,1:11)){
 
#text(10,120-(i*10),paste("cluster",i,"ctrl vs spdh :",mvec[i]),pos=4)
text(10,130-(i*10),mvecWT[i],pos=4,col=rgb(0,64,255,maxColorValue=256))
text(28,130-(i*10),rep("/",12),pos=4)
text(30,130-(i*10),mvecSPDH[i],pos=4,col=rgb(255,64,0,maxColorValue=256))
text(60,130-(i*10),mvec2[i],pos=4)
}

plot(1,xlim=c(1,100),ylim=c(1,130),col="white",main="max peak value",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
for(i in c(12,1:11)){

#text(10,120-(i*10),paste("cluster",i,"ctrl vs spdh :",mvec[i]),pos=4)
text(10,130-(i*10),mvecwt0[i],pos=4,col=rgb(0,64,255,maxColorValue=256))
text(28,130-(i*10),rep("/",12),pos=4)
text(30,130-(i*10),mvecspdh0[i],pos=4,col=rgb(255,64,0,maxColorValue=256))
text(60,130-(i*10),mvec2[i],pos=4)
}
dev.off()


## plot the numbers now
pdf(paste(gene,"_numbersonly.pdf",sep=""),12,6)
par(mfrow=c(1,2))

plot(1,xlim=c(1,100),ylim=c(1,130),col="white",main="mean values of regions",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
for(i in c(12,1:11)){
 
#text(10,120-(i*10),paste("cluster",i,"ctrl vs spdh :",mvec[i]),pos=4)
text(29,130-(i*10),mvecWT[i],pos=2,col=rgb(0,64,255,maxColorValue=256))
text(29,130-(i*10),rep("/",12),pos=4)
text(30,130-(i*10),mvecSPDH[i],pos=4,col=rgb(255,64,0,maxColorValue=256))
text(60,130-(i*10),mvec2[i],pos=4)
}

plot(1,xlim=c(1,100),ylim=c(1,130),col="white",main="max peak value",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
for(i in c(12,1:11)){

#text(10,120-(i*10),paste("cluster",i,"ctrl vs spdh :",mvec[i]),pos=4)
text(28,130-(i*10),mvecwt0[i],pos=2,col=rgb(0,64,255,maxColorValue=256))
text(29,130-(i*10),rep("/",12),pos=4)
text(30,130-(i*10),mvecspdh0[i],pos=4,col=rgb(255,64,0,maxColorValue=256))
text(60,130-(i*10),mvec2[i],pos=4)
}
dev.off()





#!/usr/bin/Rscript

## you need to  have the pbmc.ctrl object already, pbmc.spdh and pbmc.combined objects
## you can generate it by running the script first in an R session 
## Figure5/ScRNA-Seq_data_processing/labbook_analyses.R
## and then run everything following below


mycolors=c("#a65628","#c49a6c","#add8e6","#00008b","#377eb8","#ff7f00","#fbb040","#f9ed32","#4daf4a","#f881bf","#e41a1c")
avx2 <- AverageExpression(pbmc.combined)

p <- pbmc.ctrl
current.cluster.ids=0:10
new.cluster.ids=1:11
p@ident<- plyr::mapvalues(p@ident, from = current.cluster.ids, to = new.cluster.ids)


cluster.n_ctrl <- table(p@ident)
cluster.n_spdh <- table(spdh_avex.cl)


par(oma=c(2,2,1,1),mar=c(5,5,2,1),mfrow=c(3,4),cex=1.2)

t1t=p@dr$tsne@cell.embeddings

plot(t1t[,1],t1t[,2],xlab="t-SNE 1",ylab="t-SNE 2",main="",col="white",xlim=c(-40,80),las=1)
for(k in 1:11){
    points(t1t[p@ident==k,1],t1t[p@ident==k,2],col=mycolors[k],pch=16)
}
legend("topright",legend=1:11,col=mycolors,pch=16,bty="n",cex=1.2,xjust=1)


#mycolors=c("#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404", "#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")
t2=t1t
t2[,2]=(t2[,2]+40)/50
t2[,1]=(((t2[,1]+40)/80)*2)+3

for(i in 1:11){
    j=i-1
    ## markers is on scale 0:10
    names_markers <- markers[[1]][which(markers[[1]]$cluster==j),7]
    hl=which(rownames(leftavex)%in%names_markers)
    hl=which(rownames(avx2)%in%names_markers)

    plot(log(avx2[,i]+1),log(avx2[,i+11]+1),
         main=paste("Ctrl cluster",i),
         xlab=paste("Average expression [#Umis]\nCTRL -",cluster.n_ctrl[i],"cells"),
         ylab=paste("Average expression [#Umis]\nSPDH - assigned to cluster",i," - ",cluster.n_spdh[i],"cells"),
         col="grey",
         pch=16,
         xlim=c(0,10),
         ylim=c(0,10),
         xaxt="n"
         )
    abline(0,1,col="red")
    axis(1,c(0,2,4,6,8),c(0,2,4,6,8))
    points(log(avx2[hl,i]+1),log(avx2[hl,i+1]+1),col="blue",pch=16)

    points(t2[,1]+5,t2[,2],col="grey")
    points(t2[p@ident==i,1]+5,t2[p@ident==i,2],col=mycolors[i],pch=16)
}




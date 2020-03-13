#!/usr/bin/Rscript

source("dplot_function.R")
## you need to  have the pbmc.ctrl object already
## you can generate it by running the script first in an R session 
## Figure5/ScRNA-Seq_data_processing/labbook_analyses.R
## and then run everything following below

p <- pbmc.ctrl
current.cluster.ids=0:10
new.cluster.ids=1:11
p@ident<- plyr::mapvalues(p@ident, from = current.cluster.ids, to = new.cluster.ids)





genesin <-c("Pcna","Twist1","Prrx2","Dlk1","Col1a1","Hoxa11","Meox2","Shox2","Msx1","Msx2","Hoxd12","Aldh1a2","Irx1","Inhba","Sox6","Foxc1","Col9a2","Fibin","Col2a1","Matn4","Pax9","Myod1","Myog","Tnnt2","Fcgr1","Cd33","Cd48","Hbb-bh1","Hbb-bs","Hba-a2","Hist1h2bc")

df.wt <- getmat(p=p,g=genesin,clustl=11)
df.wt$Cluster=factor(df.wt$Cluster)
df.wt$Cluster2 <- factor(df.wt$Cluster, levels = c("8", "11", "10", "9", "7", "6","4","5","3","2","1"))

pl1 <- ggplot(df.wt, aes(x=Gene, y=Cluster2, color=Frac, size=Exprs)) +geom_point(alpha=0.6)+
scale_x_discrete(limit = genesin)+ labs(title = "Genes top10 markers clusters 1:8", x = "", y = "Clusters")+scale_color_gradient( high = "#132B43", low = "#56B1F7")+
theme(axis.text.x=element_text(angle = 45, hjust = +1))+
scale_size_continuous(range = c(0.1,10),breaks=seq(0,1,0.25))

save_plot(paste("FigureS4B.pdf",sep=""),plot_grid(pl1,ncol=1,nrow=1),nrow=1,ncol=1,base_height=10,base_width=20)



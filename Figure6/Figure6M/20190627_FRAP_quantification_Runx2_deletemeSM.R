#1. Load required packages and read in data
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(plyr)
library(Rmisc)
library(data.table)

setwd("~/Desktop/FRAP_data_quantification")
FRAP_Runx2 <- read.csv("FRAP_merge_Runx2-series.csv")

summary_FRAP_Runx2 <- summarySE(FRAP_Runx2, measurevar="Normalized_intensity", groupvars=c("Timepoint","Phenotype"))


#2. Required plotting function
plot_FRAP_Runx2_paper <- function(dataframe) {
  summary <- summarySE(dataframe, measurevar="Normalized_intensity", groupvars=c("Timepoint","Phenotype"))
  summary_paper <- filter(summary, summary$Phenotype != "3Q")
  summary_paper$Normalized_intensity <- as.numeric(summary_paper$Normalized_intensity)
  summary_paper$Timepoint<- as.numeric(summary_paper$Timepoint)
  plot <- ggplot(summary_paper,aes(x=Timepoint,
                                        y=Normalized_intensity,
                                        colour=Phenotype)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=Normalized_intensity-sd, ymax=Normalized_intensity+sd),width=0.2) +
    scale_color_manual(values=c("#a97c50","#662d91"),labels=c("Runx2 +10A","Runx2 wt")) +
    scale_x_discrete(limits=c("1",
                              "2",
                              "3",
                              "4",
                              "5",
                              "6",
                              "7",
                              "8",
                              "9",
                              "10",
                              "11",
                              "12",
                              "13",
                              "14",
                              "15")) +
    labs(
      x="Time post bleach (sec)",
      y="Fluorescence intensity (%)",
      color= " "
    ) +
    theme_bw(base_size=10) +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=11),
          legend.title = element_text(size=11),
          text=element_text(family="Helvetica"),
          panel.grid.major=element_line(color="grey",linetype = "dashed",size=0.1),
          panel.grid.minor=element_line(color="grey",linetype = "dashed",size=0.1)) 
  plot
}
## This is the function for Figure 6M

#3. Save FRAP plot comparison between different constructs (Figure 6M)
plot_FRAP_Runx2_paper(FRAP_Runx2)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="FRAP of the Runx2 series(wt, +10A) with mean value and standard deviation - paper version.pdf",
       plot=last_plot(),  
       width=20,
       height=12.5,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

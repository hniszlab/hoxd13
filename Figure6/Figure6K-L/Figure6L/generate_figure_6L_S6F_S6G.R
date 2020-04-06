#!/usr/bin/Rscript

# generate_figure_6L_S6F_S6G.R -  R script
# Copyright (C) 2018 Dora Knezevic
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#1. Read in data from excel files as a list with each track as an individual data frame within list
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(plyr)
library(Rmisc)
library(data.table)

colnames <- c("Min, Intensities #1",
              "Max, Intensities #1",
              "Mean, Intensities #1",
              "Sum, Intensities #1",
              "SD, Intensities #1",
              "SNR (Mean/SD), Intensities #1",
              "Min, Intensities #2",
              "Max, Intensities #2",
              "Mean, Intensities #2",
              "Sum, Intensities #2",
              "SD, Intensities #2",
              "SNR (Mean/SD), Intensities #2",
              "Area, Projection (XY/Z)",
              "Area, Projection (XY/Z) without holes")
# here the column names are extracted from an example data frame, 
# but this should be refined to be automatically extracted from the df

#1.1 Required functions
Intensity0_vs_phase_shift <- function(dataframe) {
  timepoint1 <- filter(dataframe,dataframe[,17] == 1)
  intensity_timepoint1 <- select(timepoint1, `Mean, Intensities #1`,`Cell track`)
  phaseshift_intensity1 <- merge(dataframe, intensity_timepoint1, by="Cell track", all = TRUE)
  phaseshift_intensity1$`Mean, Intensities #1.y` <- as.numeric(phaseshift_intensity1$`Mean, Intensities #1.y`)
  phaseshift_intensity1
}


phase_shift_score <- function(area_with_droplets,area_without_droplets) {
  area_with_droplets <- as.numeric(area_with_droplets)
  area_without_droplets <- as.numeric(area_without_droplets)
  ps_score <- (area_without_droplets-area_with_droplets)/area_without_droplets
  ps_score
} 

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  sheets <- sheets[-1] # this deletes the master first sheet
  x <- lapply(sheets, function(X) readxl::read_xlsx(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x <- lapply(x,setNames,colnames) # this sets appropriate names to the columns (has to be predefined)
  x <- lapply(x, function(x) x[-1,]) # this deletes the first row with the repeated column names
  x <- x[sapply(x,nrow) == 10] # this selects only the data frames which have at least 10 rows i.e measured across 10 timepoints
  x <- lapply(x, cbind, Timepoint = c(1:10)) #this annotates timepoint so it can be extracted later
  x <- bind_rows(x, .id = "Cell track") # this concatenes all the data frames (cell tracks) into one dataframe per image 
  x$"Phase-shift score" <- phase_shift_score(as.numeric(x[,14]), as.numeric(x[,15])) #this creates a phase-shift score column
  x
}
# if an image has no cells (tracks) with 10 points this will register as an empty dataframe 
# later on and will cause an error in the function  

add_phenotype_info <- function(images) {
  images$Phenotype <- substr(str_split_fixed(images[,1], "_", 5)[,4],1,10) #this extracts only the phenotype from the character of strings
  images
}

intensity_plot_Runx2 <- function(images) {
  initial <- filter(images, images[,17] == 1)
  initial_nomCh <- filter(initial, initial[,19] != "1")
  plot <- ggplot(initial_nomCh, aes(x=initial_nomCh[,19], y=as.numeric(initial_nomCh[,5]),fill=initial_nomCh[,19],color=initial_nomCh[,19])) + 
    geom_violin(trim=F, alpha=0.4) +
    scale_fill_manual(values=c("#db222a","#e89005"),labels=c("mCh-Cry2","Runx2 wt")) + 
    scale_color_manual(values=c("#db222a","#e89005")) + 
    guides(color = FALSE) +
    geom_boxplot(width=0.1, fill="white",color=c("#db222a","#e89005")) +
    labs(
      x=" ",
      y="Mean intensity at timepoint 1",
      fill=" ") +
    theme_bw(base_size=10) +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=11),
          legend.title = element_text(size=11),
          legend.position = "right",
          text=element_text(family="Helvetica"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) 
  plot
}

intensity_plot_Runx2_paper <- function(images) {
  initial <- filter(images, images[,17] == 1)
  initial_nomCh <- filter(initial, initial[,19] != "1")
  plot <- ggplot(initial_nomCh, aes(x=initial_nomCh[,19], y=as.numeric(initial_nomCh[,5]), color =initial_nomCh[,19] )) + 
    geom_dotplot(binaxis = "y", binwidth = 1, stackdir = "center", fill= "white", position = position_jitter(width = 0.0, height = 0.9)) +
    stat_summary(fun.data = "mean_sdl", colour = "black", geom = "crossbar", width = 0.5) +
    scale_y_continuous(limits=c(18,70)) +
    scale_color_manual(values=c("#808285","#662D91"),labels=c("mCh-Cry2","Runx2 wt")) + 
    labs(
      x=" ",
      y="Mean intensity at timepoint 1") +
    theme_bw(base_size=10) +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size=11),
          legend.title = element_text(size=11),
          legend.position = "right",
          text=element_text(family="Helvetica"),
          panel.grid.major=element_line(color="grey",linetype = "dashed",size=0.1),
          panel.grid.minor=element_line(color="grey",linetype = "dashed",size=0.1)) 
  plot
}
## This is the function for Figure S6F

intensity_density_plot_Runx2 <- function(images) {
  initial <- filter(images, images[,17] == 1)
  initial_nomCh <- filter(initial, initial[,19] != "1")
  plot <- ggplot(initial_nomCh, aes(x=as.numeric(initial_nomCh[,5]),color=initial_nomCh[,19],fill=initial_nomCh[,19])) + 
    geom_density(alpha=0.4) +
    geom_rug(show.legend=F) +
    #geom_histogram(aes(y=..density..), alpha=0.5, position="identity",bins = 10) +
    scale_color_manual(values=c("#db222a","#e89005")) +
    scale_fill_manual(values=c("#db222a","#e89005"),labels=c("mCh-Cry2","Runx2 wt")) + 
    guides(color = FALSE) +
    labs(
      x="Intensity at timepoint 1",
      y="Density",
      fill=" ") +
    theme_bw(base_size=10) +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=11),
          legend.title = element_text(size=11),
          text=element_text(family="Helvetica"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) 
  plot
}

plot_phase_shift_Runx2_paper <- function(dataframe) {
  dataframe_nomCh <- filter(dataframe, dataframe[,19] != "1")
  summary <- summarySE(dataframe_nomCh, measurevar="Phase-shift score", groupvars=c("Timepoint","Phenotype"))
  summary$Score <- summary$"Phase-shift score"
  plot <- ggplot(summary,aes(x=as.numeric(Timepoint),
                             y=Score,
                             colour=Phenotype)) +
    geom_point() +
    geom_errorbar(aes(ymin=Score-se, ymax=Score+se),width=0.2) +
    geom_line() +
    scale_color_manual(values=c("#808285","#662D91"),labels=c("mCh-Cry2","Runx2 wt")) +
    scale_x_discrete(limits=c("1",
                              "2",
                              "3",
                              "4",
                              "5",
                              "6",
                              "7",
                              "8",
                              "9",
                              "10")) + # this needs to be fixed - why does it not sort 1-10
    labs(
      x="Timepoint",
      y="Phase-shift score",
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
## This is the function for Figure 6L

phase_shift_score_density_plot <- function(images) {
  initial <- filter(images, images[,17] == 1)
  initial_nomCh <- filter(initial, initial[,19] != "1")
  plot <- ggplot(initial_nomCh, aes(x=initial_nomCh[,18],color=initial_nomCh[,19],fill=initial_nomCh[,19])) + 
    geom_density(alpha=0.4) +
    geom_rug(show.legend=F) +
    #geom_histogram(aes(y=..density..), alpha=0.5, position="identity",bins = 20) +
    scale_color_manual(values=c("#1F7A8C","#A31621")) +
    scale_fill_manual(values=c("#1F7A8C","#A31621"),labels=c("Hoxd13 DEdel","Hoxd13 wt")) + 
    guides(color = FALSE) +
    labs(
      x="Phase-shift score at timepoint 1",
      y="Density",
      fill=" ") +
    theme_bw(base_size=10) +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=11),
          legend.title = element_text(size=11),
          text=element_text(family="Helvetica"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) 
  plot
}

intensity_DDC_plot_Runx2 <- function(dataframe) {
  dataframe_mod <- Intensity0_vs_phase_shift(dataframe)
  dataframe_nomCh <- filter(dataframe_mod, dataframe_mod[,19] != "1")
  plot <- ggplot(dataframe_nomCh,aes(x=dataframe_nomCh[,20],
                                     y=dataframe_nomCh$"Phase-shift score",
                                     colour=Phenotype)) +
    scale_y_continuous(limits=c(0,0.2)) +
    geom_point() +
    facet_wrap(~ Timepoint, nrow=2, ncol=5) +
    scale_color_manual(values=c("#db222a","#e89005"),labels=c("mCh-Cry2","Runx2 wt")) +
    labs(
      x="Intensity",
      y="Phase-shift score",
      color= " ") +
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

intensity_DDC_plot_Runx2_only10 <- function(dataframe) {
  dataframe_mod <- Intensity0_vs_phase_shift(dataframe)
  dataframe_nomCh <- filter(dataframe_mod, dataframe_mod[,19] != "mCh-Cry2")
  dataframe_nomCh_10 <- filter(dataframe_nomCh, dataframe_nomCh[,17] == 10)
  plot <- ggplot(dataframe_nomCh_10,aes(x=dataframe_nomCh_10[,20],
                                        y=dataframe_nomCh_10$"Phase-shift score",
                                        colour=Phenotype)) +
    scale_y_continuous(limits=c(0,0.10)) +
    scale_x_continuous(limits=c(0,70)) +
    geom_point(size=4) +
    stat_cor(method = "pearson") + 
    scale_color_manual(values=c("#662D91"),labels=c("Runx2 wt")) +
    labs(
      x="Intensity",
      y="Phase-shift score",
      color= " ") +
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
## This is the function for Figure S6G 

cwd=getwd()
#1.2 Load the required files and run function over them  
setwd("../../../Fdata/Raw data from the image quantification - Runx2wt and mCh-Cry2")
Droplet_quantification <- list.files(pattern="xlsx$")

image_file <- vector("list")

for (i in Droplet_quantification){
  image_name <- substr(i,1,40) 
  print(image_name)
  image_file[[paste(image_name, sep="_")]] <- read_excel_allsheets(i)
}

all_images_Runx2 <- bind_rows(image_file, .id="df")

all_images_Runx2 <- add_phenotype_info(all_images_Runx2)

#1.3 Gating and setting cut-offs 
# 1.3.1 Removing cells which "phase separate" (FP) at timepoint 1
cells_timepoint_0_Runx2 <- filter(all_images_Runx2,all_images_Runx2[,17] == 1)

no_ps_at_t0_cells_Runx2 <- filter(cells_timepoint_0_Runx2, cells_timepoint_0_Runx2[,18] == 0)
no_ps_at_t0_cells_Runx2$`Mean, Intensities #1` <- as.numeric(no_ps_at_t0_cells_Runx2$`Mean, Intensities #1`)

all_images_filtered_Runx2 <- subset(all_images_Runx2, all_images_Runx2$`Cell track` %in% no_ps_at_t0_cells_Runx2$`Cell track`)

#1.3.2 Correcting for intensity differences between populations
fullhist = hist(as.numeric(no_ps_at_t0_cells_Runx2$`Mean, Intensities #1`), breaks = 20)
wthist = with(subset(no_ps_at_t0_cells_Runx2, no_ps_at_t0_cells_Runx2$Phenotype == "wt"), hist(as.numeric(`Mean, Intensities #1`), breaks = fullhist$breaks))
condition_hist = with(subset(no_ps_at_t0_cells_Runx2, no_ps_at_t0_cells_Runx2$Phenotype == "mCh-Cry2"), hist(as.numeric(`Mean, Intensities #1`), breaks = fullhist$breaks))

combhist = fullhist
combhist$counts = wthist$counts - condition_hist$counts
plot(combhist)

number_of_cells_to_remove_Runx2 <- combhist$counts
number_of_cells_to_remove_Runx2
# a positive number indicates the cells should be removed from wt, a negative one that they should be removed from condition

set.seed(123)

bin_20_25 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 20,25))
bin_20_25 <- data.table(bin_20_25)
bin_20_25_filtered <- bin_20_25[-sample(which(bin_20_25$Phenotype=="wt"), 29)]

bin_25_30 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 25,30))
bin_25_30 <- data.table(bin_25_30)
bin_25_30_filtered <- bin_25_30[-sample(which(bin_25_30$Phenotype=="wt"), 37)]

bin_30_35 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 30,35))
bin_30_35 <- data.table(bin_30_35)
bin_30_35_filtered <- bin_30_35[-sample(which(bin_30_35$Phenotype=="wt"), 68)]

bin_35_40 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 35,40))
bin_35_40 <- data.table(bin_35_40)
bin_35_40_filtered <- bin_35_40[-sample(which(bin_35_40$Phenotype=="wt"), 14)]

bin_40_45 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 40,45))
bin_40_45 <- data.table(bin_40_45)
bin_40_45_filtered <- bin_40_45#[-sample(which(bin_40_45$Phenotype=="wt"), 0)]

bin_45_50 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 45,50))
bin_45_50 <- data.table(bin_45_50)
bin_45_50_filtered <- bin_45_50[-sample(which(bin_45_50$Phenotype=="mCh-Cry2"), 19)]

bin_50_55 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 50,55))
bin_50_55 <- data.table(bin_50_55)
bin_50_55_filtered <- bin_50_55[-sample(which(bin_50_55$Phenotype=="mCh-Cry2"), 26)]

bin_55_60 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 55,60))
bin_55_60 <- data.table(bin_55_60)
bin_55_60_filtered <- bin_55_60[-sample(which(bin_55_60$Phenotype=="mCh-Cry2"), 11)]

bin_60_65 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 60,65))
bin_60_65 <- data.table(bin_60_65)
bin_60_65_filtered <- bin_60_65[-sample(which(bin_60_65$Phenotype=="mCh-Cry2"), 13)]

bin_65_70 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 65,70))
bin_65_70 <- data.table(bin_65_70)
bin_65_70_filtered <- bin_65_70[-sample(which(bin_65_70$Phenotype=="mCh-Cry2"), 16)]

bin_70_75 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 70,75))
bin_70_75 <- data.table(bin_70_75)
bin_70_75_filtered <- bin_70_75[-sample(which(bin_70_75$Phenotype=="mCh-Cry2"), 14)]

bin_75_80 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 75,80))
bin_75_80 <- data.table(bin_75_80)
bin_75_80_filtered <- bin_75_80[-sample(which(bin_75_80$Phenotype=="mCh-Cry2"), 15)]

bin_80_85 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 80,85))
bin_80_85 <- data.table(bin_80_85)
bin_80_85_filtered <- bin_80_85[-sample(which(bin_80_85$Phenotype=="mCh-Cry2"), 10)]

bin_85_90 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 85,90))
bin_85_90 <- data.table(bin_85_90)
bin_85_90_filtered <- bin_85_90[-sample(which(bin_85_90$Phenotype=="mCh-Cry2"), 13)]

bin_90_95 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 90,95))
bin_90_95 <- data.table(bin_90_95)
bin_90_95_filtered <- bin_90_95[-sample(which(bin_90_95$Phenotype=="mCh-Cry2"), 7)]

bin_95_100 <- filter(no_ps_at_t0_cells_Runx2,between(no_ps_at_t0_cells_Runx2[,5], 95,100))
bin_95_100 <- data.table(bin_95_100)
bin_95_100_filtered <- bin_95_100[-sample(which(bin_95_100$Phenotype=="mCh-Cry2"), 3)]


intensity_filtered_Runx2 <- rbind(bin_20_25_filtered,
                                  bin_25_30_filtered,
                                  bin_30_35_filtered,
                                  bin_35_40_filtered,
                                  bin_40_45_filtered,
                                bin_45_50_filtered,
                                bin_50_55_filtered,
                                bin_55_60_filtered,
                                bin_60_65_filtered,
                                bin_65_70_filtered,
                                bin_70_75_filtered,
                                bin_75_80_filtered,
                                bin_80_85_filtered,
                                bin_85_90_filtered,
                                bin_90_95_filtered,
                                bin_95_100_filtered)

all_images_filtered_intensity_Runx2 <- subset(all_images_Runx2, all_images_Runx2$`Cell track` %in% intensity_filtered_Runx2$`Cell track`)


setwd(cwd)
# 2. Plot relevant data
# 2.0 Intensity dot plot paper style (Figure S6F)
intensity_plot_Runx2_paper(all_images_filtered_intensity_Runx2)
ggsave(filename="FigureS6F.pdf",
       plot=last_plot(),   
       width=17,
       height=13,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 2.3 Phase-shift plot with mean and standard error of the mean (Figure 6L)
plot_phase_shift_Runx2_paper(all_images_filtered_intensity_Runx2)

ggsave(filename="Figure6L.pdf",
       plot=last_plot(),   
       width=20,
       height=12.5,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 2.6 Dose-dependent curve intensity vs. phase-shift score in only phase 10 (Figure S6G)
intensity_DDC_plot_Runx2_only10(all_images_filtered_intensity_Runx2)
ggsave(filename="FigureS6G.pdf",
       plot=last_plot(),  
       width=20,
       height=15,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

#!/usr/bin/Rscript

# helper_functions.R -  R script
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
# this is specific to the way I name my files, but can be adapted by changing ,4 to other number depending on word ordering

intensity_plot <- function(images) {
  initial <- filter(images, images[,17] == 1) # takes only timepoint 1 measurements
  initial_nomCh <- filter(initial, initial[,19] != "1") # removes mCherry (labelled as 1 in my phenotype column)
  plot <- ggplot(initial_nomCh, aes(x=initial_nomCh[,19], y=as.numeric(initial_nomCh[,5]),fill=initial_nomCh[,19],color=initial_nomCh[,19])) + 
    geom_violin(trim=F, alpha=0.4) +
    scale_fill_manual(values=c("#1F7A8C","#A31621"),labels=c("Hoxd13 DEdel","Hoxd13 wt")) + 
    scale_color_manual(values=c("#1F7A8C","#A31621")) + 
    guides(color = FALSE) +
    geom_boxplot(width=0.1, fill="white",color=c("#1F7A8C","#A31621")) +
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

intensity_density_plot <- function(images) {
  initial <- filter(images, images[,17] == 1)
  initial_nomCh <- filter(initial, initial[,19] != "1")
  plot <- ggplot(initial_nomCh, aes(x=as.numeric(initial_nomCh[,5]),color=initial_nomCh[,19],fill=initial_nomCh[,19])) + 
    geom_density(alpha=0.4) +
    geom_rug(show.legend=F) +
    #geom_histogram(aes(y=..density..), alpha=0.5, position="identity",bins = 10) +
    scale_color_manual(values=c("#1F7A8C","#A31621")) +
    scale_fill_manual(values=c("#1F7A8C","#A31621"),labels=c("Hoxd13 DEdel","Hoxd13 wt")) + 
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

plot_phase_shift <- function(dataframe) {
  dataframe_nomCh <- filter(dataframe, dataframe[,19] != "1") 
  summary <- summarySE(dataframe_nomCh, measurevar="Phase-shift score", groupvars=c("Timepoint","Phenotype"))
  summary$Score <- summary$"Phase-shift score"
  plot <- ggplot(summary,aes(x=as.numeric(Timepoint),
                               y=Score,
                               colour=Phenotype)) +
    geom_point() +
    geom_errorbar(aes(ymin=Score-se, ymax=Score+se),width=0.2) +
    geom_smooth(method="loess", se=F) +
    scale_color_manual(values=c("#1F7A8C","#A31621"),labels=c("Hoxd13 DEdel","Hoxd13 wt")) +
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

ecdf_plot <- function(dataframe){
  dataframe_nomCh <- filter(dataframe, dataframe[,19] != "1")
  plot <- ggplot(dataframe_nomCh, aes(as.numeric(dataframe_nomCh[,18]),colour=dataframe_nomCh[,19])) + 
    stat_ecdf(geom = "step") +
    scale_color_manual(values=c("#1F7A8C","#A31621"),labels=c("Hoxd13 DEdel","Hoxd13 wt")) +
    labs(
      x="Phase-shift score",
      y="F(x)",
      colour = " "
    ) +
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

Intensity0_vs_phase_shift <- function(dataframe) {
  timepoint1 <- filter(dataframe,dataframe[,17] == 1)
  intensity_timepoint1 <- select(timepoint1, `Mean, Intensities #1`,`Cell track`)
  phaseshift_intensity1 <- merge(dataframe, intensity_timepoint1, by="Cell track", all = TRUE)
  phaseshift_intensity1$`Mean, Intensities #1.y` <- as.numeric(phaseshift_intensity1$`Mean, Intensities #1.y`)
  phaseshift_intensity1
}

intensity_DDC_plot <- function(dataframe) {
  dataframe_mod <- Intensity0_vs_phase_shift(dataframe) #this uses the function from above
  dataframe_nomCh <- filter(dataframe_mod, dataframe_mod[,19] != "1")
  plot <- ggplot(dataframe_nomCh,aes(x=dataframe_nomCh[,20],
                             y=dataframe_nomCh$"Phase-shift score",
                             colour=Phenotype)) +
    scale_y_continuous(limits=c(0,0.2)) +
    geom_point() +
    facet_wrap(~ Timepoint, nrow=2, ncol=5) +
    scale_color_manual(values=c("#1F7A8C","#A31621"),labels=c("Hoxd13 DEdel","Hoxd13 wt")) +
    labs(
      x="Intensity at 0''",
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

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

#1.2 Load the required files and run function over them  
setwd("~/Desktop/Master thesis/Droplet quantification/DEdel/2-4") ##folder name is changed depending on the condition read in
Droplet_quantification <- list.files()

image_file <- vector("list")

for (i in Droplet_quantification){
  image_name <- substr(i,1,40) 
  print(image_name)
  image_file[[paste(image_name, sep="_")]] <- read_excel_allsheets(i)
}

all_images <- bind_rows(image_file, .id="df")

all_images <- add_phenotype_info(all_images)

#1.3 Gating and setting cut-offs 
# 1.3.1 Removing cells which "phase separate" (FP) at timepoint 1
cells_timepoint_0 <- filter(all_images,all_images[,17] == 1)

no_ps_at_t0_cells <- filter(cells_timepoint_0, cells_timepoint_0[,18] == 0)
no_ps_at_t0_cells$`Mean, Intensities #1` <- as.numeric(no_ps_at_t0_cells$`Mean, Intensities #1`)


all_images_filtered <- subset(all_images, all_images$`Cell track` %in% no_ps_at_t0_cells$`Cell track`)

#1.3.2 Correcting for intensity differences between populations
fullhist = hist(as.numeric(no_ps_at_t0_cells$`Mean, Intensities #1`), breaks = 20)
wthist = with(subset(no_ps_at_t0_cells, no_ps_at_t0_cells$Phenotype == "wt"), hist(as.numeric(`Mean, Intensities #1`), breaks = fullhist$breaks))
condition_hist = with(subset(no_ps_at_t0_cells, no_ps_at_t0_cells$Phenotype == "DEdel"), hist(as.numeric(`Mean, Intensities #1`), breaks = fullhist$breaks))

combhist = fullhist
combhist$counts = wthist$counts - condition_hist$counts
plot(combhist)

number_of_cells_to_remove <- combhist$counts
number_of_cells_to_remove
# a positive number indicates the cells should be removed from wt, a negative one that they should be removed from condition

## these are numbers removed from DEdel condition, the numbers are adjusted manually depending on what condition is read in
set.seed(123)
bin_30_35 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 30,35))
bin_30_35 <- data.table(bin_30_35)
bin_30_35_filtered <- bin_30_35[-sample(which(bin_30_35$Phenotype=="wt"), 1)]

bin_35_40 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 35,40))
bin_35_40 <- data.table(bin_35_40)
bin_35_40_filtered <- bin_35_40[-sample(which(bin_35_40$Phenotype=="DEdel"), 3)]

bin_40_45 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 40,45))
bin_40_45 <- data.table(bin_40_45)
bin_40_45_filtered <- bin_40_45[-sample(which(bin_40_45$Phenotype=="wt"), 3)]

bin_45_50 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 45,50))
bin_45_50 <- data.table(bin_45_50)
bin_45_50_filtered <- bin_45_50[-sample(which(bin_45_50$Phenotype=="DEdel"), 4)]

bin_50_55 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 50,55))
bin_50_55 <- data.table(bin_50_55)
bin_50_55_filtered <- bin_50_55[-sample(which(bin_50_55$Phenotype=="wt"), 1)]

bin_55_60 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 55,60))
bin_55_60 <- data.table(bin_55_60)
bin_55_60_filtered <- bin_55_60[-sample(which(bin_55_60$Phenotype=="DEdel"), 4)]

bin_60_65 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 60,65))
bin_60_65 <- data.table(bin_60_65)
bin_60_65_filtered <- bin_60_65[-sample(which(bin_60_65$Phenotype=="DEdel"), 4)]

bin_65_70 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 65,70))
bin_65_70 <- data.table(bin_65_70)
bin_65_70_filtered <- bin_65_70[-sample(which(bin_65_70$Phenotype=="wt"), 2)]

bin_70_75 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 70,75))
bin_70_75 <- data.table(bin_70_75)
bin_70_75_filtered <- bin_70_75[-sample(which(bin_70_75$Phenotype=="wt"), 5)]

bin_75_80 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 75,80))
bin_75_80 <- data.table(bin_75_80)
bin_75_80_filtered <- bin_75_80[-sample(which(bin_75_80$Phenotype=="wt"), 6)]

bin_80_85 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 80,85))
bin_80_85 <- data.table(bin_80_85)
bin_80_85_filtered <- bin_80_85#[-sample(which(bin_80_85$Phenotype=="wt"), 0)]

bin_85_90 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 85,90))
bin_85_90 <- data.table(bin_85_90)
bin_85_90_filtered <- bin_85_90[-sample(which(bin_85_90$Phenotype=="DEdel"), 4)]

bin_90_95 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 90,95))
bin_90_95 <- data.table(bin_90_95)
bin_90_95_filtered <- bin_90_95[-sample(which(bin_90_95$Phenotype=="wt"), 3)]

bin_95_100 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 95,100))
bin_95_100 <- data.table(bin_95_100)
bin_95_100_filtered <- bin_95_100[-sample(which(bin_95_100$Phenotype=="wt"), 5)]

bin_100_105 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 100,105))
bin_100_105 <- data.table(bin_100_105)
bin_100_105_filtered <- bin_100_105#[-sample(which(bin_100_105$Phenotype=="-15A"), 5)]

bin_105_110 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 105,110))
bin_105_110 <- data.table(bin_105_110)
bin_105_110_filtered <- bin_105_110[-sample(which(bin_105_110$Phenotype=="wt"), 4)]

bin_110_115 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 110,115))
bin_110_115 <- data.table(bin_110_115)
bin_110_115_filtered <- bin_110_115[-sample(which(bin_110_115$Phenotype=="wt"), 3)]

bin_115_120 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 115,120))
bin_115_120 <- data.table(bin_115_120)
bin_115_120_filtered <- bin_115_120[-sample(which(bin_115_120$Phenotype=="DEdel"), 1)]

bin_120_125 <- filter(no_ps_at_t0_cells,between(no_ps_at_t0_cells[,5], 120,125))
bin_120_125 <- data.table(bin_120_125)
bin_120_125_filtered <- bin_120_125[-sample(which(bin_120_125$Phenotype=="wt"), 1)]

intensity_filtered <- rbind(bin_30_35_filtered,
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
                                bin_95_100_filtered,
                                bin_100_105_filtered,
                                bin_105_110_filtered,
                                bin_110_115_filtered,
                                bin_115_120_filtered,
                                bin_120_125_filtered)

all_images_filtered_intensity <- subset(all_images, all_images$`Cell track` %in% intensity_filtered$`Cell track`)

# 2. Plot relevant data  
# 2.1 Violin plot of intensity distribution across conditions
intensity_plot(all_images_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Intensity violin plot of Hoxd13wt vs. Hoxd13DEdel after correction.pdf",
       plot=last_plot(),   
       width=17,
       height=13,
       units = "cm",
       dpi=600)

dev.off()

# 2.2 Density distribution curve of intensity across conditions
intensity_density_plot(all_images_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Density distribution plot of Hoxd13wt vs. Hoxd13DEdel after correction.pdf",
       plot=last_plot(),   
       width=17,
       height=13,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 2.3 Phase-shift plot with mean and standard error of the mean
plot_phase_shift(all_images_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Phase-shift between conditions mean + standard error of the mean Hoxd13wt vs. Hoxd13DEdel after intensity correction with fitted curve.pdf",
       plot=last_plot(),   
       width=20,
       height=12.5,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 2.4 ECDF plot of the phase-shift score
ecdf_plot(all_images_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="ECDF plot of the phase-shift score in Hoxd13wt vs. Hoxd13DEdel after intensity correction.pdf",
       plot=last_plot(),
       width=17,
       height=13,
       units = "cm",
       dpi=600)

dev.off()


# 2.5 Dose-dependent curve intensity vs. phase-shift score 
intensity_DDC_plot(all_images_filtered_intensity) 

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Concentration(intensity) dependent effect on the phase-shift score in Hoxd13wt vs. Hoxd13DEdel after intensity correction.pdf",
       plot=last_plot(),  
       width=20,
       height=10,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()


# 3. Summary statistics overview
summary_DEdel <- summarySE(all_images_filtered_intensity, measurevar="Phase-shift score", groupvars=c("Timepoint","Phenotype"))

droplet_forming_cells_DEdel <- filter(all_images_filtered_intensity, all_images_filtered_intensity[,18] > 0)
summary_droplet_forming_cells_DEdel <- summarySE(droplet_forming_cells_DEdel, measurevar="Phase-shift score", groupvars=c("Timepoint","Phenotype"))


all_images_filtered_intensity$`Mean, Intensities #1` <- as.numeric(all_images_filtered_intensity$`Mean, Intensities #1`)

# 3.1 Wilcox test and KS test of mean intensities between conditions post correction 
wilcox.test(all_images_filtered_intensity[all_images_filtered_intensity[,19] == "DEdel",5], all_images_filtered_intensity[all_images_filtered_intensity[,19] == "wt",5], paired = FALSE, alternative = "two.sided")

ks.test(all_images_filtered_intensity[all_images_filtered_intensity[,19] == "DEdel",5], all_images_filtered_intensity[all_images_filtered_intensity[,19] == "wt",5])


# 4. Saving the output
setwd("~/Desktop/Master thesis/Quantification_final") 
write.csv(all_images_filtered_intensity, "201906_Hoxd13DEdel-plates-analysis_series2-4_DK.csv")
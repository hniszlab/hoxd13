##Note: Subsections starting with X are not relevant for the final figure generation.

# 0. Necessary functions
plot_phase_shift_all_constructs <- function(dataframe) {
  summary <- summarySE(dataframe, measurevar="Phase-shift score", groupvars=c("Timepoint","Phenotype"))
  summary$Score <- summary$"Phase-shift score"
  plot <- ggplot(summary,aes(x=as.numeric(Timepoint),
                             y=Score,
                             colour=Phenotype)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Score-se, ymax=Score+se),width=0.2) +
    scale_color_manual(values=c("#19647E", "#119DA4","#009E5E","#662D91"),labels=c("Hoxd13 -15A","Hoxd13 -7A","Hoxd13 DEdel" ,"Hoxd13 wt")) +
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
## This is the function for 7F

intensity_plot_paper <- function(images) {
  initial <- filter(images, images[,17] == 1)
  plot <- ggplot(initial, aes(x=initial[,19], y=as.numeric(initial[,5]), color =initial[,19] )) + 
    geom_dotplot(binaxis = "y", binwidth = 1, stackdir = "center", fill= "white", position = position_jitter(width = 0.0, height = 0.9)) +
    stat_summary(fun.data = "mean_sdl", colour = "black", geom = "crossbar", width = 0.5) +
    scale_y_continuous(limits=c(70,115)) +
    scale_color_manual(values=c("#f6924e", "#e6424a","#4777be","#981349"),labels=c("Hoxd13 -15A","Hoxd13 -7A","Hoxd13 DEdel" ,"Hoxd13 wt")) +
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
## This is the function for S7F

Parameter_dependency_plot <- function(dataframe) {
  initial <- filter(dataframe, dataframe$Timepoint == 10) 
  plot <- ggplot(initial,aes(x=Ala,
                             y=`Phase-shift score`,
                             colour=Phenotype)) +
    scale_y_continuous(limits=c(0,0.1)) +
    scale_x_continuous(limits=c(24,36)) +
    geom_point(size=4) +
    geom_errorbar(aes(ymin=`Phase-shift score`-se, ymax=`Phase-shift score`+se),width=0.4) +
    stat_cor(method = "spearman") + 
    scale_color_manual(values=c("#19647E", "#119DA4","#009E5E","#662D91"),labels=c("Hoxd13 -15A","Hoxd13 -7A","Hoxd13 DEdel" ,"Hoxd13 wt")) +
    labs(
      x="% Ala",
      y="Phase-shift score at 180''",
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
# this is only for timepoint 10 phase-shift score
## This is the function for 7G

intensity_DDC_plot_all_constructs_only10 <- function(dataframe) {
  dataframe_mod <- Intensity0_vs_phase_shift(dataframe)
  dataframe_mod_10 <- filter(dataframe_mod, dataframe_mod[,17] == 10)
  plot <- ggplot(dataframe_mod_10,aes(x=dataframe_mod_10[,20],
                                      y=dataframe_mod_10$"Phase-shift score",
                                      colour=Phenotype)) +
    scale_y_continuous(limits=c(0,0.15)) +
    scale_x_continuous(limits=c(70,105)) +
    geom_point(size=4) +
    stat_cor(method = "pearson") + 
    scale_color_manual(values=c("#f6924e", "#e6424a","#4777be","#981349"),labels=c("Hoxd13 -15A","Hoxd13 -7A","Hoxd13 DEdel" ,"Hoxd13 wt")) +
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
## This is the function for S7G

# 1. Read in data (here it was already read in so only a duplicate was made)
all_images_15A_2 <- all_images_15A
all_images_7A_2 <- all_images_7A
all_images_DEdel_2 <- all_images
## The origin of this data/function to extract it is in the 190516_Droplet_quantification_ArivisVision4D.R
## section 1., 1.1. and 1.2. (change working directory folder depending on condition, necessary raw data is included)

# X. Re-annotate wt so it is associated to its parent test condition
all_images_15A_2$Phenotype[(all_images_15A_2$Phenotype) == "wt"] <- "wt_15A"
all_images_7A_2$Phenotype[(all_images_7A_2$Phenotype) == "wt"] <- "wt_7A"
all_images_DEdel_2$Phenotype[(all_images_DEdel_2$Phenotype) == "wt"] <- "wt_DEdel"
##This was performed only to check intervariation between wt conditions, but wt was not differentiated
#  in the final figure. 


# 2. Combine all conditions
all_constructs <- rbind(all_images_15A_2,
                        all_images_7A_2,
                        all_images_DEdel_2)


# 3. Select only cells in a specific intensity range to ensure comparability 

# 3.1. Filter all cells which have a phase-shift score >0 at timepoint 1
cells_timepoint_0_all <- filter(all_constructs,all_constructs[,17] == 1)
cells_timepoint_0_all$`Mean, Intensities #1` <- as.numeric(cells_timepoint_0_all$`Mean, Intensities #1`)

no_ps_at_t0_cells_all <- filter(cells_timepoint_0_all, cells_timepoint_0_all[,18] == 0)
no_ps_at_t0_cells_all$`Mean, Intensities #1` <- as.numeric(no_ps_at_t0_cells_all$`Mean, Intensities #1`)

all_constructs_filtered <- subset(all_constructs, all_constructs$`Cell track` %in% no_ps_at_t0_cells_all$`Cell track`)

# 3.2 Select 30 cells from each condition with intensity range of timepoint 1 between 75 and 100
bin_75_100 <- filter(no_ps_at_t0_cells_all,between(no_ps_at_t0_cells_all[,5], 75,100))
bin_75_100 <- data.table(bin_75_100)
bin_75_100_wt <- bin_75_100[sample(which(bin_75_100$Phenotype=="wt"), 30)] ## here wt is banded together from all three conditions
bin_75_100_DEdel <- bin_75_100[sample(which(bin_75_100$Phenotype=="DEdel"), 30)]
bin_75_100_7A <- bin_75_100[sample(which(bin_75_100$Phenotype=="-7A"), 30)]
bin_75_100_15A <- bin_75_100[sample(which(bin_75_100$Phenotype=="-15A"), 30)]

intensity_filtered_all <- rbind(bin_75_100_wt,
                                bin_75_100_DEdel,
                                bin_75_100_7A,
                                bin_75_100_15A)

# X.2 Select 30 cells from each condition with intensity range of timepoint 1 between 75 and 100 without combining wt
bin_75_100 <- filter(no_ps_at_t0_cells_all,between(no_ps_at_t0_cells_all[,5], 75,100))
bin_75_100 <- data.table(bin_75_100)
bin_75_100_wt_DEdel <- bin_75_100[sample(which(bin_75_100$Phenotype=="wt_DEdel"), 30)] # here wt is separated across conditions so inter-wt differences can be observed
bin_75_100_wt_7A <- bin_75_100[sample(which(bin_75_100$Phenotype=="wt_7A"), 30)]
bin_75_100_wt_15A <- bin_75_100[sample(which(bin_75_100$Phenotype=="wt_15A"), 30)]
bin_75_100_DEdel <- bin_75_100[sample(which(bin_75_100$Phenotype=="DEdel"), 30)]
bin_75_100_7A <- bin_75_100[sample(which(bin_75_100$Phenotype=="-7A"), 30)]
bin_75_100_15A <- bin_75_100[sample(which(bin_75_100$Phenotype=="-15A"), 30)]

intensity_filtered_all_full <- rbind(bin_75_100_wt_DEdel,
                                     bin_75_100_wt_7A,
                                     bin_75_100_wt_15A,
                                     bin_75_100_DEdel,
                                     bin_75_100_7A,
                                     bin_75_100_15A)

# 3.3 Merge full information from relevant 30 cells
all_constructs_filtered_intensity <- subset(all_constructs, all_constructs$`Cell track` %in% intensity_filtered_all$`Cell track`)

all_constructs_filtered_intensity_full <- subset(all_constructs, all_constructs$`Cell track` %in% intensity_filtered_all$`Cell track`)

# 3.4 Check if intensity distribution is significantly different between conditions
wilcox.test(as.numeric(bin_75_100_15A$`Mean, Intensities #1`),as.numeric(bin_75_100_DEdel$`Mean, Intensities #1`))
wilcox.test(as.numeric(bin_75_100_15A$`Mean, Intensities #1`),as.numeric(bin_75_100_7A$`Mean, Intensities #1`))
wilcox.test(as.numeric(bin_75_100_15A$`Mean, Intensities #1`),as.numeric(bin_75_100_wt$`Mean, Intensities #1`))
wilcox.test(as.numeric(bin_75_100_7A$`Mean, Intensities #1`),as.numeric(bin_75_100_DEdel$`Mean, Intensities #1`))
wilcox.test(as.numeric(bin_75_100_7A$`Mean, Intensities #1`),as.numeric(bin_75_100_wt$`Mean, Intensities #1`))
wilcox.test(as.numeric(bin_75_100_DEdel$`Mean, Intensities #1`),as.numeric(bin_75_100_wt$`Mean, Intensities #1`))


# 4. Create dataframe summarizing additional construct characteristics together with mean phase shift score 
## This is relevant for Figure 7G
summary <- summarySE(all_constructs_filtered_intensity, measurevar="Phase-shift score", groupvars=c("Timepoint","Phenotype"))
summary$GRAVY <- 1
summary$Ala <- 1
summary$DE <- 1
summary$pI <- 1

summary$GRAVY[(summary$Phenotype) == "wt"] <- "0.090"
summary$GRAVY[(summary$Phenotype) == "-15A"] <- "-0.145"
summary$GRAVY[(summary$Phenotype) == "-7A"] <- "-0.012"
summary$GRAVY[(summary$Phenotype) == "DEdel"] <- "0.273"

summary$Ala[(summary$Phenotype) == "wt"] <- "33.9"
summary$Ala[(summary$Phenotype) == "-15A"] <- "24.8"
summary$Ala[(summary$Phenotype) == "-7A"] <- "29.9"
summary$Ala[(summary$Phenotype) == "DEdel"] <- "35.6"

summary$DE[(summary$Phenotype) == "wt"] <- "4.8"
summary$DE[(summary$Phenotype) == "-15A"] <- "5.6"
summary$DE[(summary$Phenotype) == "-7A"] <- "5.2"
summary$DE[(summary$Phenotype) == "DEdel"] <- "0"

summary$pI[(summary$Phenotype) == "wt"] <- "7.84"
summary$pI[(summary$Phenotype) == "-15A"] <- "7.84"
summary$pI[(summary$Phenotype) == "-7A"] <- "7.84"
summary$pI[(summary$Phenotype) == "DEdel"] <- "11.3"

summary$`Phase-shift score` <- as.numeric(summary$`Phase-shift score`)
summary$GRAVY <- as.numeric(summary$GRAVY)
summary$Ala <- as.numeric(summary$Ala)
summary$DE <- as.numeric(summary$DE)
summary$pI <- as.numeric(summary$pI)


# 6. Plot graphical representations 
# 6.1. Dose-dependent curve parameters vs. phase-shift score (Figure 7G)
Parameter_dependency_plot(summary)
## which parameter is plotted should be changed in the function by specifying different column

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Dependency of phase-shift score on  percent Ala in the HOXD13 deletion and substitution series.pdf",
       plot=last_plot(),  
       width=17,
       height=13,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 6.2. Intensity at timepoint 1 plot (Figure S7F)
intensity_plot_paper(all_constructs_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Intensity dotplot plot of Hoxd13 wt vs. -7A vs- -15A vs. DEdel after intensity gating for 30 cells between 75 and 100 with mean and standard deviation.pdf",
       plot=last_plot(),   
       width=17,
       height=13,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 6.3. Phase-shift score across time all construct combined (Figure 7F)
plot_phase_shift_all_constructs(all_constructs_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Phase-shift between conditions mean + standard error of the mean Hoxd13 wt vs. -7A vs. -15A vs. DEdel after intensity gating for 30 cells between 75 and 100 with connected points curve - paper version.pdf",
       plot=last_plot(),   
       width=20,
       height=12.5,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

# 6.4. Dose-dependent curve intensity vs. phase-shift score in only phase 10 (Figure S7G)
intensity_DDC_plot_all_constructs_only10(all_constructs_filtered_intensity)

setwd("~/Desktop/Master thesis/Figures") 
ggsave(filename="Concentration(intensity) dependent effect on the phase-shift score in the Hoxd13 series after intensity correction only at timepoint 10 with Pearson correlation.pdf",
       plot=last_plot(),  
       width=25,
       height=15,
       units = "cm",
       dpi=600,
       useDingbats=FALSE)

dev.off()

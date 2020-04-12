# Copyright (C) 2018  Shaon Basu
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

#START OF Dual Flouresecence Plot (DFP)

#setwd
setwd()
#specify path to folder containing analysis files


#read data
red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("mix2", red[i,5])==TRUE)
  {red[i,6]="d13(7A) mix 2"}
  else if(grepl("d13wt", red[i,5])==TRUE)
  {red[i,6]="d13(WT)"}
  else if(grepl("2mM", red[i,5])==TRUE)
  {red[i,6]="d13(7A) + \n2 mM ATP"}
  else if(grepl("-4mM", red[i,5])==TRUE)
  {red[i,6]="d13(7A) + \n4 mM ATP"}
  else if(grepl("8mM", red[i,5])==TRUE)
  {red[i,6]="d13(7A) + \n8 mM ATP"}
  else if(grepl("16mM", red[i,5])==TRUE)
  {red[i,6]="d13(7A) + \n16 mM ATP"}
  else if(grepl("24mM", red[i,5])==TRUE)
  {red[i,6]="d13(7A) + \n24 mM ATP"}
  else{red[i,6]="d13(7A)"}
}
#View(red)
names(red)[6] <- "mix"
#green$group <- gsub("-.*","",green$file)
red$mix<- factor(red$mix, levels = c("d13(7A)", "d13(7A) + \n2 mM ATP","d13(7A) + \n4 mM ATP","d13(7A) + \n8 mM ATP","d13(7A) + \n16 mM ATP","d13(7A) + \n24 mM ATP", "d13(WT)", "d13(7A) mix 2"))
red2 <- subset(red, mix == "d13(7A)" | mix == "d13(7A) + \n16 mM ATP" | mix == "d13(7A) + \n8 mM ATP" | mix == "d13(WT)")
red1 <- subset(red, mix == "d13(7A)" | mix == "d13(7A) + \n24 mM ATP" | mix == "d13(7A) + \n16 mM ATP" | mix == "d13(7A) + \n8 mM ATP" | mix == "d13(7A) + \n4 mM ATP" | mix == "d13(7A) + \n2 mM ATP" | mix == "d13(WT)")

#make plots
require(ggplot2)
require(RColorBrewer)

#make Fig3I
ggplot(red1, aes(x=mix, y=green/red)) + theme_classic() + geom_jitter() + ylab('Med1 content / D13 content') + scale_y_log10()

ggplot(red1, aes(x=mix, y=green/red)) + theme_classic() + geom_boxplot(outlier.shape = NA) + ylab('Med1 content / D13 content') + scale_y_log10()

#make Fig#H
ggplot(red2, aes(y=red, x=green, color=mix)) + geom_point(alpha=0.4, size = red2$area/500000) + theme_classic() + ylab('HOXD13 condensates [mCherry signal]') + xlab('Mediator enrichment [GFP signal]') + scale_color_manual(values=c('red','darkorange','lightseagreen','blue')) + scale_x_log10()


#END OF DFP


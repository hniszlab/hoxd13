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

#setwd to Figure 6Q Image analysis folder
setwd()

#read values

red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("Runx2_0,75uM", red[i,5])==TRUE)
  {red[i,6]="Runx2WT 0,75 uM"}
  else if(grepl("Runx2_1,5uM", red[i,5])==TRUE)
  {red[i,6]="Runx2WT 1,5 uM"}
  else if(grepl("10A)_1,5uM_Med1_3uM", red[i,5])==TRUE)
  {red[i,6]="Runx2(10A) 1,5 uM"}
  else{red[i,6]="Runx2(10A) 0,75 uM"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix, levels = c("Runx2WT 0,75 uM", "Runx2WT 1,5 uM", "Runx2(10A) 1,5 uM","Runx2(10A) 0,75 uM"))


runx2 <- subset(red, mix == "Runx2WT 0,75 uM" | mix == "Runx2(10A) 0,75 uM")


#make Figure 6Q
ggplot(runx2, aes(y=red, x=green, color=mix)) + geom_point(alpha=0.4, size = runx2$area/500000) + theme_classic() + ylab('HOXD13 condensates [mCherry signal]') + xlab('Mediator enrichment [GFP signal]') + scale_color_manual(values=c('blue','red'))



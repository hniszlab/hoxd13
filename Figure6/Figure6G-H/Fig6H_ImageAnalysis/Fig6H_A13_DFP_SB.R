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


#setwd to folder containing 6H Image analysis files
setwd("")

#read table of values
red <- read.csv2("objectsRED", header=T, sep = ",")
red$file <- as.character(red$file)
for(i in 1: nrow(red))
{
  if(grepl("A13_3,5uM", red[i,5])==TRUE)
  {red[i,6]="Hoxa13WT 3,5 uM"}
  else{red[i,6]="Hoxa13(+7A) 3,5uM"}
}
names(red)[6] <- "mix"
red$mix <- factor(red$mix)

#make Figure 6H Dual Flourescence Plot
ggplot(red, aes(y=red, x=green, color=mix)) + geom_point(alpha=0.4, size = red$area/500000) + theme_classic() + ylab('HOXD13 condensates [mCherry signal]') + xlab('Mediator enrichment [GFP signal]') + scale_color_manual(values=c('red','blue'))






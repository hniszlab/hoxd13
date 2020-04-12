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



#set wd, to folder "Fig1C_VutaraSRX_ClusterAnalysis" 
setwd('')

#read data
gatedcell <- read.csv2("Fig1C-TableofValues.csv", header=T, sep = ",", dec = ".")

cell2 <- subset(gatedcell, gatedcell$ID == 'cell2')

#Fig 1C Inset
require(ggplot2)
require(RColorBrewer)

ggplot(cell2, aes(x=size, color = ID, fill = ID)) + theme_classic() + xlab("2 x Radius of Gyration (nm)") +
  geom_density(alpha=.4) + scale_fill_manual(values=c('purple','blue')) + scale_color_manual(values=c('purple','blue')) +
  ggtitle("STORM : Hoxd13 Condensate Size  \nParticle Count 15 - 50")

summary(cell2$size)



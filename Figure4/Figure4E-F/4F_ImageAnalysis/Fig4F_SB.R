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


#setwd, specify path to folder containing analysis files ('4F_ImageAnalysis')

setwd()


#read values

merge3 <- read.csv("Hoxd13tether_LacIO_TotalMerge_Gated_200115.csv", sep = ",")
merge4 <- subset(merge3, merge3$mix == "Med1" | merge3$mix == "Hoxd13 WT" | merge3$mix == "Hoxd13 7A" | merge3$mix == "Hoxd13 10A")


#Make Fig 4F

ggplot(merge4, aes(x=mix, y=(merge4$TetherYellow/merge4$RingYellow))) + theme_classic() + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) + ylab('yellow enrichment at tether') + ylim(0.9,1.6) 


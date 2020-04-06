# generate_figure_S7H.R -  R script
# Copyright (C) 2019 Henri Niskanen
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

library(ggplot2)

# Read in data
df <- read.delim("FigureS7H_raw_data.txt", sep = "\t", header=TRUE, stringsAsFactors = FALSE)

# exclude GAL4-DBD control values before plotting
df2 <- subset(df, sample!="GAL4-DBD")
df2$sample_as_factor <- as.factor(df2$sample)

# Plot values
p <- ggplot(df2, aes(x=sample_as_factor, y=log_FC)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Print plot
pdf("Figure_S7H.pdf")
plot(p, width = 10, height = 20)
dev.off()
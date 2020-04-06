#!/usr/bin/Rscript

# generate_figure_6R.R -  R script
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
library(car)

# Read in data
df <- read.delim("Figure6R_raw_data.txt", sep = "\t", header=TRUE, stringsAsFactors = FALSE)

# exclude GAL4-DBD control values before plotting
df2 <- subset(df, sample!="GAL4-DBD")
df2$sample_as_factor <- as.factor(df2$sample)

# Plot values
p <- ggplot(df2, aes(x=sample_as_factor, y=log_FC)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Values to test
a <- df2[df2$sample=="RUNX2wt","log_FC"]
b <- df2[df2$sample=="RUNX2(+10A)","log_FC"]

# Test normality
shapiro.test(a)
# Shapiro-Wilk normality test
# data:  a
# W = 0.99012, p-value = 0.9895
shapiro.test(b)
# Shapiro-Wilk normality test
# data:  b
# W = 0.90621, p-value = 0.4119

# Test variance
group <- as.factor(c(rep(1, length(a)), rep(2, length(b))))
leveneTest(c(a,b), group)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value  Pr(>F)  
# group  1  6.2125 0.03185 *
#       10    

# Welch Two Sample t-test
t.test(a,b)
# 
# data:  a and b
# t = 25.458, df = 7.0634, p-value = 3.279e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   3.067072 3.693904
# sample estimates:
#   mean of x mean of y 
# 1.940586 -1.439902 

# Print plot
pdf("Figure_6R.pdf", width = 3, height = 10)
plot(p)
dev.off()

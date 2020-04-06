#!/usr/bin/Rscript

# generate_figure_7H-I.R -  R script
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

#install.packages("ggplot2")
library(ggplot2)

#https://github.com/kassambara/easyGgplot2
#install.packages("devtools")
#library(devtools)
#install_github("kassambara/easyGgplot2")
library(easyGgplot2)
library(car)

###############################################
# Plot Figure 7H
###############################################

# Read in data
df <- read.delim("../../Fdata/Figure7H-I/Figure7H_raw_data.txt", sep = "\t", header=TRUE, stringsAsFactors = FALSE)

# exclude GAL4-DBD control values before plotting
df2 <- subset(df, sample!="GAL4-DBD")
df2$sample_as_factor <- as.factor(df2$sample)


# MNX1 sample comparisons

# Values to test
a <- df2[df2$sample=="MNX1wt","log_FC"]
b <- df2[df2$sample=="MNX1-Adel","log_FC"]


# Test normality
shapiro.test(a)
# Shapiro-Wilk normality test
# data:  a
# W = 0.89975, p-value = 0.3725
shapiro.test(b)
# Shapiro-Wilk normality test
# data:  b
# W = 0.888, p-value = 0.1903

# Test variance
group <- as.factor(c(rep(1, length(a)), rep(2, length(b))))
leveneTest(c(a,b), group)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value  Pr(>F)  
# group  1  5.5836 0.03438 *
#       13                  

# Welch Two Sample t-test
t.test(a,b)
# data:  a and b
# t = -3.341, df = 11.583, p-value = 0.006152
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.8076906 -0.1685086
# sample estimates:
#   mean of x mean of y 
# 2.425152  2.913251 


# OLIG2 sample comparisons

# Values to test
a <- df2[df2$sample=="OLIG2wt","log_FC"]
b <- df2[df2$sample=="OLIG2-Adel","log_FC"]

# Test normality
shapiro.test(a)
# Shapiro-Wilk normality test
# data:  a
# W = 0.89665, p-value = 0.3545
shapiro.test(b)
# Shapiro-Wilk normality test
# data:  b
# W = 0.94889, p-value = 0.6779

# Test variance
group <- as.factor(c(rep(1, length(a)), rep(2, length(b))))
leveneTest(c(a,b), group)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value  Pr(>F)  
# group  1  3.5296 0.08289 .
#       13                  

# Two Sample t-test
t.test(a,b, var.equal = TRUE)
# data:  a and b
# t = -7.4453, df = 13, p-value = 4.866e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.184231 -1.201744
# sample estimates:
#   mean of x mean of y 
# 1.165086  2.858073 


t.test(a,b)
# Welch Two Sample t-test
# data:  a and b
# t = -8.6462, df = 11.54, p-value = 2.218e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.121508 -1.264467
# sample estimates:
#   mean of x mean of y 
# 1.165086  2.858073 



# TBX1 sample comparisons

# Values to test
a <- df2[df2$sample=="TBX1wt","log_FC"]
b <- df2[df2$sample=="TBX1-Adel","log_FC"]

# Test normality
shapiro.test(a)
# Shapiro-Wilk normality test
# data:  a
# W = 0.91047, p-value = 0.1377
shapiro.test(b)
# Shapiro-Wilk normality test
# data:  b
# W = 0.91625, p-value = 0.1688

# Test variance
group <- as.factor(c(rep(1, length(a)), rep(2, length(b))))
leveneTest(c(a,b), group)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  1   0.011 0.9173
#       28                

# Two Sample t-test
t.test(a,b, var.equal = TRUE)
# data:  a and b
# t = -4.6263, df = 28, p-value = 7.701e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.7871879 -0.6902373
# sample estimates:
#   mean of x mean of y 
# 2.110050  3.348762 


t.test(a,b)
# Welch Two Sample t-test
# data:  a and b
# t = -4.6263, df = 27.999, p-value = 7.702e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.7871887 -0.6902365
# sample estimates:
#   mean of x mean of y 
# 2.110050  3.348762 




# EOMES sample comparisons

# Values to test
a <- df2[df2$sample=="EOMESwt","log_FC"]
b <- df2[df2$sample=="EOMES-Adel","log_FC"]

# Test normality
shapiro.test(a)
# Shapiro-Wilk normality test
# data:  a
# W = 0.85807, p-value = 0.04625
# --> not normally distributed

shapiro.test(b)
# Shapiro-Wilk normality test
# data:  b
# W = 0.78181, p-value = 0.005863
# --> not normally distributed


# Test variance
group <- as.factor(c(rep(1, length(a)), rep(2, length(b))))
leveneTest(c(a,b), group)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  1  0.1712  0.683
#       22              

# Two Sample t-test
t.test(a,b, var.equal = TRUE)
# Two Sample t-test
# data:  a and b
# t = -0.63688, df = 22, p-value = 0.5308
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.5596303  0.2966645
# sample estimates:
#   mean of x mean of y 
# 6.366797  6.498280 


t.test(a,b)

# --> p-value > 0.05 in EOMESwt vs EOMES-Adel,
# also, data not normally distributed






# HOXA13 sample comparisons

# Values to test
a <- df2[df2$sample=="HOXA13wt","log_FC"]
b <- df2[df2$sample=="HOXA13-Adel","log_FC"]

# Test normality
shapiro.test(a)
# Shapiro-Wilk normality test
# data:  a
# W = 0.9627, p-value = 0.7393
shapiro.test(b)
# Shapiro-Wilk normality test
# data:  b
# W = 0.94218, p-value = 0.4106

# Test variance
group <- as.factor(c(rep(1, length(a)), rep(2, length(b))))
leveneTest(c(a,b), group)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value   Pr(>F)   
# group  1  10.474 0.003105 **
#       28                        

# Welch Two sample t-test
t.test(a,b)
# data:  a and b
# t = -7.7313, df = 20.477, p-value = 1.68e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.4141978 -0.8139276
# sample estimates:
#   mean of x mean of y 
# 0.2290013 1.3430640 

# Plot values
p <- ggplot(df2, aes(x=sample_as_factor, y=log_FC)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Print plot
pdf("Figure_7H.pdf")
plot(p, width = 10, height = 20)
dev.off()




###############################################
# Plot Figure 7I
###############################################


df <- read.delim("../../Fdata/Figure7H-I/Figure7I_raw_data.txt", sep = "\t", header=TRUE, stringsAsFactors = FALSE)

df2 <- subset(df, sample!="GAL4-DBD")

plot(df2$A_content, y=df2$log_FC)

cor.test(df2$A_content, df2$log_FC, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  df2$A_content and df2$log_FC
# t = -14.213, df = 73, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9074984 -0.7822698
# sample estimates:
#        cor 
# -0.8570639 

abline(lm(log_FC ~ A_content, data=df2))  # plot regline
linearMod <- lm(log_FC ~ A_content, data=df2)  # build linear regression model on full data
summary(linearMod)

# Call:
#   lm(formula = log_FC ~ A_content, data = df2)
#
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.8685 -0.2624 -0.0327  0.3343  0.7422 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.631620   0.156973   16.77   <2e-16 ***
#   A_content   -0.064278   0.004522  -14.21   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3929 on 73 degrees of freedom
# Multiple R-squared:  0.7346,	Adjusted R-squared:  0.7309 
# F-statistic:   202 on 1 and 73 DF,  p-value: < 2.2e-16


# Plot
p <- ggplot2.scatterplot(data=df2, xName='A_content',yName='log_FC', 
                    #groupName='sample', 
                    size=2,
                    backgroundColor="white", removePanelGrid = TRUE, 
                    addRegLine = TRUE, addConfidenceInterval = TRUE)


# Print plot
pdf("Figure_7I.pdf")
plot(p, width = 20, height = 20)
dev.off()




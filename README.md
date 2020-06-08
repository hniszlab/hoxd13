# hoxd13 `Readme`

## About 
### Unblending of Transcriptional Condensates in Human Repeat Expansion Disease. 
#### Published:May 07, 2020 DOI:https://doi.org/10.1016/j.cell.2020.04.018

Authors: Sebastian Mackowiak, Shaon Basu, Henri Niskanen, Dora Knezevic

This github contains mainly the computational analysis to the study. 

Complementary data for running most of the R scripts can be found at Mendeley, doi: 10.17632/ztd6wzcv7h.13

## Getting a copy 

1a. You can checkout this repository by typing on the command line

```
git clone https://github.com/hniszlab/hoxd13.git
```

## or

1b. download the zip file and extract it by typing
``` 
wget https://github.com/hniszlab/hoxd13/archive/master.zip -O hoxd13.zip
unzip hoxd13.zip
```

2. Change the directory of the hoxd13 project
```
cd hoxd13
```


## Running the code
All R code and packages were run in R version 3.6.0 (2019-04-26)

Start an R console and install the necessary dependencies if they should be missing
```
R
library(dendextend)
library(cowplot)
library(svglite)
library(gplots)
library(Matrix)
library(Rtsne)
library(irlba)
library(raster)
library(RColorBrewer)
library(scales)
library(data.table)
library(argparse)
library(directlabels)
```
Missing R libraries can be installed by the install.packages command.
E.g.
```
install.packages("Matrix")
```

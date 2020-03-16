# hoxd13 `Readme`

## About 

Author: Sebastian Mackowiak 
The scripts found in this repository were used to generate
various figures in ... 

Complementary data for the scripts can be found at
mendeley.com/hniszlab/Fdata

## Getting a copy 

1a. First checkout this github by typing on the command line

```
git clone git@github.com/hnisz/hoxd13
```

## or

1b. download the zip file and extract it by typing
``` 
unzip hoxd13.zipo
```

2. Change the directory of the hoxd13 project
```
cd hodx13
```

3. Get the complementary data by typing
```
wget mendeley.com/hnisz/Fdata
```

## Running the code
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
```
Missing R libraries can be installed by the install.packages command.
E.g.
```
install.packages("Matrix")
```

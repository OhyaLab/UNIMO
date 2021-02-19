# Calculating the best smooth span value (f) for LOESS regression using
# lo() function of 'gamlss' package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46)
# and Akaike information criterion (AIC).

## Of the 501 CalMorph parameters (Ohya et al, 2005; PNAS 102(52): 19015-19020), 220 parameters are
## the coefficient of variation (CV) of their related mean values. CV parameters are non-linearly
## dependent on the mean values. To uncouple this concomitant dependency, LOESS (locally
## estimated scatterplot smoothing) regression with a smooth span (f) has been reported
## effective (Levy and Siegal, 2008; PLoS Biol. 6(11): e264).
## Finding the best value of 'f' can be challenging. 'f' can be set between 0.10 to 0.99.
## Here, we employed gamlss() function to test various 'f' values: f = {0.10, 0.11, 0.12, 0.13, ..., 0.99}.
## Eventually, the best 'f' value is chosen by AIC.

# Authors: Ghanegolmohammadi et al. (2021) 
##################################################
##################################################
##################################################
## Setting the working directory (not necessary to define).
setwd("...")
## Setting the initial seed (not necessary to set).
set.seed(123)

# Prerequisite
## 'gamlss' package is used for generalized linear modeling.
### Looking for the 'gamlss' package among the installed ones.
if("gamlss" %in% rownames(installed.packages())){
  ### Loading 'gamlss' package into R environment.
  library(gamlss)
  ### If 'gamlss' package is not installed, installing and loading it.
} else {install.packages("gamlss"); library(gamlss)}

# Main object to import the data of 114 diploid Wild-type yeast strains.
Diploid_WT <- NULL
## Loading 114 diploid Wild-type yeast strains from 'SCMD2: Saccharomyces Cerevisiae Morphological Database 2':
### Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114data.tsv
Diploid_WT$data <- read.table("wt114data.tsv", header = TRUE, row.names = 1)
#### While importing the data,  R sometimes converts '-' to '.' in column names.
#### To have consistency in parameter names, this problem should be fixed
colnames(Diploid_WT$data) <- chartr(".", "-", colnames(Diploid_WT$data))

## Defining a vector of mean parameters with corresponding 'CV' parameter as its names.
idx <- colnames(Diploid_WT$data)[grep("CV", colnames(Diploid_WT$data))]
names(idx) <- idx
idx <- sub("CV", "", idx)

# A matrix to save the best 'f' value and its AIC.
## This matrix has 220 CV parameters as its row-names.
SmoothSpan <- matrix(data = NA, nrow = length(idx), ncol = 2, dimnames = list(names(idx), c("Smooth span", "AIC")))

for(i in names(idx)){
  temp <- NULL
  aics <- NULL
  ## Defining a data frame with mean and CV parameters:
  ### x = Mean parameters.
  ### y = CV parameters.
  temp <- data.frame(x = Diploid_WT$data[,idx[i]], y = Diploid_WT$data[,i])
  ## Fitting a LOESS regression given various 'f' values: f = {0.10, 0.11, 0.12, 0.13, ..., 0.99}
  ## and checking goodness of fit using AIC() function of stats package
  ## LOESS regression function in 'gamlss' package is lo()
  aics <- sapply(10:99/100, function(i) AIC(gamlss(y~lo(~x, span = i), data = na.omit(temp), family = NO)))
  ## Finding the best fit using AIC --> minimum value
  SmoothSpan[i,"Smooth span"] <- (10:99/100)[which.min(aics)]
  SmoothSpan[i,"AIC"] <- min(aics)
}
## Saving the output as a R data
save(SmoothSpan, file = "SmoothSpan_WT114.rdata")

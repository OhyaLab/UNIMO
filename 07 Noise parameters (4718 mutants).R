# To test the effectiveness of the proposed pipeline, we reanalyzed morphological variations of
# the 4718 yeast nonessential gene mutants by comparing them to a dataset of haploid wild-type yeast
# strains (his3; 109 replicates) in the manner of the generalize linear model. 

## Calculating the best smooth span value (f) for LOESS regression using
## lo() function of 'gamlss' package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46)
## and Akaike information criterion (AIC).

### Of the 501 CalMorph parameters (Ohya et al, 2005; PNAS 102(52): 19015-19020), 220 parameters are
### the coefficient of variation (CV) of their related mean values. CV parameters are non-linearly
### dependent on the mean values. To uncouple this concomitant dependency, LOESS (locally
### estimated scatterplot smoothing) regression with a smooth span (f) has been reported
### effective (Levy and Siegal, 2008; PLoS Biol. 6(11): e264).
### Finding the best value of 'f' can be challenging. 'f' can be set between 0.10 to 0.99.
### Here, we employed gamlss() function to test various 'f' values: f = {0.10, 0.11, 0.12, 0.13, ..., 0.99}.
### Eventually, the best 'f' value is chosen by AIC.

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

# Main object to import data of 4718 non-essential budding yeast mutants and wild-type. 
AllHap <- NULL
# Loading 4718 non-essential budding yeast mutants:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718data.tsv
AllHap$Mut <- read.table("mt4718data.tsv", header = TRUE, row.names = 1)
### While importing the data,  R sometimes converts '-' to '.' in column names.
### To have consistency in parameter names, this problem should be fixed
colnames(AllNonEss$Mut) <- chartr(".","-", colnames(AllNonEss$Mut))

# Loading 122 haploid Wild-type (WT) strains:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122data.tsv
AllHap$WT <- read.table("wt122data.tsv", header = TRUE, row.names = 1)
### Following replicates of the WT data were not used (Suzuki et al., 2018, BMC Genomics, 19(1):149)
res <- c("04his3-12", "04his3-13", "04his3-18", "04his3-22", "04his3-33", "04his3-54", 
         "04his3-55", "04his3-60", "04his3-70", "04his385", "his3cnt-1", "his3old-1", "his3old-5")
### Deleting these replicates from the data
AllHap$WT <- AllHap$WT[!rownames(AllHap$WT) %in% res,]
rm(res)
### Final data: Aggregating WT and mutants data
AllHap$Data <- rbind(AllHap$WT, AllHap$Mut)

## Defining a vector of mean parameters with corresponding 'CV' parameter as its names.
idx <- colnames(AllHap$Data)[grep("CV", colnames(AllHap$Data))]
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
  temp <- data.frame(x=AllHap$Data[,idx[i]], y=AllHap$Data[,i])
  ## Fitting a LOESS regression given various 'f' values: f = {0.10, 0.11, 0.12, 0.13, ..., 0.99}
  ## and checking goodness of fit using AIC() function of stats package
  ## LOESS regression function in 'gamlss' package is lo()
  aics <- sapply(10:99/100, function(i) AIC(gamlss(y~lo(~x, span=i), data=temp, family=NO)))
  ## Finding the best fit using AIC --> minimum value
  SmoothSpan[i,"Smooth span"] <- (10:99/100)[which.min(aics)]
  SmoothSpan[i,"AIC"] <- min(aics)
}
rm(i)
## Saving the output as a R data
save(SmoothSpan, file = "SmoothSpan_NonEssential.rdata")
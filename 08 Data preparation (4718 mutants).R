# To test the effectiveness of the proposed pipeline, we reanalyzed morphological variations of
# the 4718 yeast nonessential gene mutants by comparing them to a data set of haploid wild-type yeast
# strains (his3; 109 replicates) in the manner of the generalize linear model.

# Transformation of the 501 CalMorph parameters (Ohya et al, 2005; PNAS 102(52): 19015-19020) including:

# 1) Uncoupling the dependency between the coefficient of variation (CV) and the mean
#    values (Levy and Siegal, 2008; PLoS Biol. 6(11): e264; Yvert et al. 2013; BMC Syst. Biol., 7(1): 54).
## Note: the best smooth span value (f) is already determined using
## lo() function of 'gamlss' package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46)
## and Akaike information criterion (AIC). See '01 Noise parameters (114 WT).r'

# 2) Handling inversion of seven ratio parameters including:
## 2-1) Invers (1/logit):
### 2-1-1) C115_A: Whole cell axis ratio at G1 phase.
### 2-1-2) C115_A1B: Mother axis ratio at S/G2 phase.
### 2-1-3) C115_C: Mother axis ratio at M phase.
## 2-2) Inverse-double (2/logit): D182_A, D182_C, D183_C, and D184_A1B
### 2-2-1) D182_A: Nuclear axis ratio at G1 phase
### 2-2-1) D182_C: Nuclear axis ratio in mother at M phase
### 2-2-1) D183_C: Nuclear axis ratio in bud at M phase
### 2-2-1) D184_A1B: Nuclear axis ratio at S/G2 phase

# Detailed information can be found in CalMorph manual at
## http://www.yeast.ib.k.u-tokyo.ac.jp/CalMorph/download.php?path=CalMorph-manual.pdf

## Authors: Ghangegolmohammadi et al. (2021)
##################################################
##################################################
##################################################

## Setting the working directory (not necessary to define).
setwd("...")
## Setting the initial seed (not necessary to set).
sset.seed(123)

# Prerequisite
## Loading an indoor LOESS regression function.
source("loess.fit.r")

## Reading information of 501 CalMorph morphological parameters.
### This file is provided along with the R scripts.
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

## Main object to save the morphological data
AllNonEss <- NULL

# Loading 4718 non-essential mutants:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718data.tsv
AllNonEss$Mut <- read.table("mt4718data.tsv", header = TRUE, row.names = 1)
colnames(AllNonEss$Mut) <- chartr(".","-", colnames(AllNonEss$Mut))
## Number of cells for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718nmrt.tsv
AllNonEss$Mut$n <- read.table("mt4718nmrt.tsv", header = TRUE, row.names = 1)
## Number of cells in specimen for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718dmnt.tsv
AllNonEss$Mut$N <- read.table("mt4718dmnt.tsv", header = TRUE, row.names = 1)

# Loading 122 haploid Wild-type strains:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122data.tsv
AllNonEss$WT <- read.table("wt122data.tsv", header = TRUE, row.names = 1)
colnames(AllNonEss$WT) <- chartr(".","-", colnames(AllNonEss$WT))
## Number of cells for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122nmrt.tsv
AllNonEss$WT$n <- read.table("wt122nmrt.tsv", header = TRUE, row.names = 1)
## Number of cells in specimen for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122dmnt.tsv
AllNonEss$WT$N <- read.table("wt122dmnt.tsv", header = TRUE, row.names = 1)

## Following replicates of the WT data were not used (Suzuki et al., 2018, BMC Genomics, 19(1):149)
res <- c("04his3-12", "04his3-13", "04his3-18", "04his3-22", "04his3-33", "04his3-54", 
         "04his3-55", "04his3-60", "04his3-70", "04his385", "his3cnt-1", "his3old-1", "his3old-5")
## Deleting these replicates from the data
AllNonEss$WT <- AllNonEss$WT[!rownames(AllNonEss$WT) %in% res,]
AllNonEss$WT$n <- AllNonEss$WT$n[!rownames(AllNonEss$WT$n) %in% res,]
AllNonEss$WT$N <- AllNonEss$WT$N[!rownames(AllNonEss$WT$N) %in% res,]
rm(res)

### Final data: Aggregating WT and mutants data
AllNonEss$Data <- rbind(AllNonEss$WT, AllNonEss$Mut)
AllNonEss$n <- rbind(AllNonEss$WT$n, AllNonEss$Mut$n)
AllNonEss$N <- rbind(AllNonEss$WT$N, AllNonEss$Mut$N)
#############################################
# Data preparation
## An object to save the transformed values
AllNonEss$Trans <- AllNonEss$Data

## Defining a vector of mean parameters with corresponding 'CV' parameter as its names.
idx <- rownames(D501)[grep("CV", rownames(D501))]
names(idx) <- idx
idx <- sub("CV", "", idx)

## Reading smooth span (defined by the lowest AIC).
### For detailed information see '01 NoiseParameters.r' file.
load("SmoothSpan_NonEssential.rdata")
## Transformtation
for(i in names(idx)) {
  ### x: Input vector of mean values.
  ### y: Input vector of CV values.
  ### ret: Return type. "all" option returns fit and noise values.
  ### f: Smooth span, to be defined in previous step (for detailed information see '01 NoiseParameters.r' script).
  temp <- lowess.fit(x = AllNonEss$Data[,idx[i]],
                     y = AllNonEss$Data[,i],
                     ret = "all",
                     f = SmoothSpan[i,"Smooth span"])
  # Noise = Observed value - Predicted value
  # Replacing the values of CV parameters by noise values.
  AllNonEss$Trans[,i] <- temp$res
}
rm(i,temp)

## Handling seven ratio parameters whose link functions are invers ("1/logit") or inver-double ("2/logit")
for(i in rownames(D501)) {
  if(D501[i,"Link_function"] == "1/logit") {
    # Inverse (reciprocal)
    AllNonEss$Trans[,i] <- 1/AllNonEss$Trans[,i]
  } else if(D501[i,"Link_function"] == "2/logit") {
    # Inverse-double (double reciprocal)
    AllNonEss$Trans[,i] <- 2/AllNonEss$Trans[,i]
  }
}
rm(i)

# Saving transformed data as a R data file.
save(AllNonEss, file = "NonEssential.Rdata")

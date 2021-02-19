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
Diploid_WT <- NULL

## Loading 114 diploid Wild-type strains:
### Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114data.tsv
Diploid_WT$Data <- read.table("wt114data.tsv", header = TRUE, row.names = 1)
colnames(Diploid_WT$Data) <- chartr(".", "-", colnames(Diploid_WT$Data))
### Number of cells for proportion parameters: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114nmrt.tsv
Diploid_WT$n <- read.table("wt114nmrt.tsv", header = TRUE, row.names = 1)
### Number of cells in specimen for proportion parameters: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114dmnt.tsv
Diploid_WT$N <- read.table("wt114dmnt.tsv", header = TRUE, row.names = 1)
#############################################
# Data preparation
## An object to save the transformed values.
Diploid_WT$Trans <- Diploid_WT$Data

## Defining a vector of mean parameters with corresponding 'CV' parameter as its names.
idx <- rownames(D501)[grep("CV", rownames(D501))]
names(idx) <- idx
idx <- sub("CV", "", idx)

## Reading smooth span (defined by the lowest AIC).
### For detailed information see '01 NoiseParameters.r' file.
load("SmoothSpan_WT114.rdata")
## Transformation
for(i in names(idx)) {
  ### x: Input vector of mean values.
  ### y: Input vector of CV values.
  ### ret: Return type. "all" option returns fit and noise values.
  ### f: Smooth span, to be defined in previous step (for detailed information see '01 NoiseParameters.r' script).
  temp <- loess.fit(x = Diploid_WT$Data[,idx[i]],
                    y = Diploid_WT$Data[,i],
                    ret = "all",
                    f = SmoothSpan[i,"Smooth span"])
  # Noise = Observed value - Predicted value
  # Replacing the values of CV parameters by noise values.
  Diploid_WT$Trans[,i] <- temp$res
}
rm(i,temp)

## Handling seven ratio parameters whose link functions are invers ("1/logit") or inver-double ("2/logit")
for(i in rownames(D501)) { 
  if(D501[i,"Link_function"] == "1/logit") {
    # Inverse (reciprocal)
    Diploid_WT$Trans[,i] <- 1/Diploid_WT$Trans[,i]
  } else if(D501[i,"Link_function"] == "2/logit") {
    # Inverse-double (double reciprocal)
    Diploid_WT$Trans[,i] <- 2/Diploid_WT$Trans[,i]
  }
}
rm(i)

# Loading confounding factor information obtained from Ohnuki and Ohya (2018; PLoS Biol. 16(5): e2005130):
## MS1: Microscope #1 before replacement of the fluorescence filter.
## MS2: Microscope #2 before replacement of the fluorescence filter.
## MS2a: Microscope #2 after replacement of the fluorescence filter over time.
## MS2b: Microscope #2 after replacement of the fluorescence filter.
## MS3: Microscope #3.
# This information will be used in modality check step.
load("Confounding_factors.Rdata")
Diploid_WT$info <- Confounding_factors

# Saving transformed data and confounding factor information as a R data file.
save(Diploid_WT, file = "Diploid_WT.Rdata")

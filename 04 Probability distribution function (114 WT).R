# Models of the probability distributions for each of the 501 CalMoprh parameters (Ohya et al, 2005; PNAS 102(52): 19015-19020)
# are needed to be determined to accommodate the statistical model used in the generalized linear modeling (GLM).
# CalMorph generates a variety of measurements that can be categorized into four groups. For each 
# group (i.e., data type), there is a well-known distribution that conventional has been used by statisticians.
# These common distributions are called reference distributions in this study.
# Rather than reference distribution for each data type, we tested various probability distributions
# to survey changes of skewness, kurtosis, and specific cases where the data are a mixture of discrete
# values (only "0" or "1", or both) and continuous values. For details see Supporting text of the paper.

## 1) Non-negative continuous values [0, +Inf) includes 183 CalMorph parameters:
### 1-1) Reference distribution: Gamma (GA)
### 1-2) Alternative distributions:
#### 1-2-1) Exponential distribution (EXP)
#### 1-2-2) Generalized Pareto distribution (GP)
#### 1-2-3) Inverse gamma distribution (IGAMMA)
#### 1-2-4) Inverse Gaussian distribution (IG)
#### 1-2-5) Log normal distribution (LOGNO)
#### 1-2-6) Log-normal 2 distribution (LOGNO2)
#### 1-2-7) Pareto type 1 distribution for y > 0 (PARETO1)
#### 1-2-8) Pareto type 2 distribution for y > 0 (PARETO2)
#### 1-2-9) Pareto type 2 distribution for y > 0 (PARETO2o)
#### 1-2-10) Weibull distribution (WEI)
#### 1-2-11) Weibull distribution for PH parameterization (WEI2)
#### 1-2-12) Weibull distribution where mu is the mean of the distribution (WEI3)
#### 1-2-13) Zero-adjusted Gamma distribution (ZAGA)
#### 1-2-14) Zero-adjusted inverse Gaussian distribution (ZAIG)

## 2) Bounded continuous values; i.e., ratio values [0,1] includes 37 CalMorph parameters:
### 2-1) Reference distribution: Beta (BE) distribution
### 2-2) Alternative distributions:
#### 2-2-1) Beta zero-inflated distribution (BEINF0)
#### 2-2-2) Beta one-inflated distribution (BEINF1)
#### 2-2-3) Beta zero- and one-inflated distribution (BEINF)
#### 2-2-4) One-inflated beta distribution (BEOI)
#### 2-2-5) Zero-inflated beta distribution (BEZI)
#### 2-2-6) Logit-normal distribution (LOGITNO)

## 3) Real continuous values; i.e., noise values (-Inf, +Inf) includes 220 CalMorph parameters:
### 3-1) Reference distribution: Gaussian (NO) distribution
### 3-2) Alternative distributions:
#### 3-2-1) Gumbel distribution (GU)
#### 3-2-2) Logistic distribution (LO)
#### 3-2-3) Normal distribution 2 with variance as sigma parameter (NO2)
#### 3-2-4) Reverse Gumbel distribution (RG)

## 4) Finite discrete values; i.e., proportion values [0,1) includes 61 CalMorph parameters:
### 4-1) Reference distribution: Binomial (BI) distribution
### 4-2) Alternative distributions:
#### 4-2-1) Beta-binomial distribution (BB)
#### 4-2-2) Zero-adjusted beta-binomial distribution (ZABB)
#### 4-2-3) Zero-adjusted binomial distribution (ZABI)
#### 4-2-4) Zero-inflated beta-binomial distribution (ZIBB)
#### 4-2-5) Zero-inflated binomial distribution (ZIBI)


# The complex mathematics of different distributions are translated into simple codes by
# 'gamlss' package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46).

# Authors: Ghanegolmohammadi et al. (2021)
##################################################
##################################################
##################################################
## Setting the working directory (not necessary to define)
setwd("...")
## Setting the initial seed (not necessary to set)
set.seed(123)

# Prerequisite
## 'gamlss' package is used for generalized linear modeling
### Looking for the 'gamlss' package among the installed ones
if("gamlss" %in% rownames(installed.packages())){
  ### Loading 'gamlss' package into R environment 
  library(gamlss)
  ### If 'gamlss' package is not installed, installing and loading it.
} else {install.packages("gamlss"); library(gamlss)}

## Loading 114 diploid Wild-type strains
### This data is already transformed and prepared in '02 Data perparation (114 WT).r' file.
load("Diploid_WT.Rdata")

## Information of the CalMorph 501 morphological parameters.
## This information is needed to find the type of the data under 'DataType' column.
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

## Controlling 'gamlss' fitting:
### Arguments are explained in gamlss.control() function of 'gamlss' package [help(gamlss.control, package = "gamlss")]:
#### n.cyc: The number of cycles of the algorithm.
#### trance: Logical; whether to print at each iteration.
gcont <- gamlss.control(n.cyc = 200L, trace = TRUE)

#########################
#Non-negative parameters#
#########################
# Finding parameters of non-negative type.
ParNames <- rownames(D501)[D501$DataType == "Non-negative"]
## Chosen probability distribution functions (PDFs) developed for continuous data ranging between [0, +Inf).
Aval_Dist <- c("EXP","GA","GP","IG","IGAMMA","LOGNO","LOGNO2","PARETO1", "PARETO2","PARETO2o",
               "WEI","WEI2","WEI3","ZAGA","ZAIG")

## An object to save fitting results.
fit <- list()

for (i in ParNames) {
  for (j in Aval_Dist) {
    temp <- NULL
    ### A data frame of non-negative values.
    temp <- data.frame(y = Diploid_WT$Trans[,i])
    ### Fitting various distributions developed for continuous data ranging between [0, +Inf).
    ### tryCatch() function outputs 'NA' in case of any error.
    fit[[i]][[j]] <-  tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = j, control = gcont), error = function(e) list(NA))
  }
}
save(fit, file = "fit(non-negative).Rdata")

## Extracting AIC values from fitted models.
### tryCatch() function outputs 'NA' in case of any error.
All_AIC <- t(sapply(ParNames, function(x) sapply(Aval_Dist, function(z) tryCatch(AIC(fit[[x]][[z]]), error = function(e) NA))))
## Finding the PDF with the minimum AIC.
Min_AIC <- apply(All_AIC, MARGIN = 1, function(x) names(which.min(x)))
## Saving the AIC and the best fitted model as a '.csv' file.
write.csv(cbind(All_AIC, Min_AIC), file = "AIC(non-negative).csv")

rm(i,j,ParNames,Aval_Dist,fit,All_AIC,Min_AIC)

#########################
#####Ratio parameters####
#########################
# Finding parameters of ratio type.
ParNames <- rownames(D501)[D501$DataType == "Ratio"]
## Chosen probability distribution functions (PDFs) developed for continuous data ranging between [0,1].
Aval_Dist <- c("BE","BEINF0","BEINF1","BEINF","BEOI","BEZI","LOGITNO","ZAGA")

### An object to save fitting results.
fit <- list()

for (i in ParNames) {
  for (j in Aval_Dist) {
    temp <- NULL
    ### A data frame of ratio values.
    temp <- data.frame(y = Diploid_WT$Trans[,i])
    ### Fitting various distributions developed for continuous data ranging between [0,1].
    ### tryCatch() function outputs 'NA' in case of any error.
    fit[[i]][[j]] <-  tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = j, control = gcont), error = function(e) list(NA))
  }
}
save(fit, file = "fit(ratio).Rdata")

## Extracting AIC values from fitted models.
### tryCatch() function outputs 'NA' in case of any error.
All_AIC <- t(sapply(ParNames, function(x) sapply(Aval_Dist, function(z) tryCatch(AIC(fit[[x]][[z]]), error = function(e) NA))))
## Finding the PDF with the minimum AIC
Min_AIC <- apply(All_AIC, MARGIN = 1, function(x) names(which.min(x)))
## Saving the AIC and the best fitted model as a '.csv' file.
write.csv(cbind(All_AIC, Min_AIC), file = "AIC(ratio).csv")

rm(i,j,ParNames,Aval_Dist,fit,All_AIC,Min_AIC)

#########################
#####Noise parameters####
#########################
# Finding parameters of noise type.
ParNames <- rownames(D501)[D501$DataType == "Noise"]
## Chosen probability distribution functions (PDFs) developed for continuous data ranging between [-Inf,+Inf].
Aval_Dist <- c("GU","LO","NO","NO2","RG")

### An object to save fitting results.
fit <- list()

for (i in ParNames) {
  for (j in Aval_Dist) {
    temp <- NULL
    ### A data frame of noise values.
    temp <- data.frame(y = Diploid_WT$Trans[,i])
    ### Fitting various distributions developed for continuous data ranging between [-Inf,+Inf].
    ### tryCatch() function outputs 'NA' in case of any error.
    fit[[i]][[j]] <-  tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = j, control = gcont), error = function(e) list(NA))
  }
}
save(fit, file = "fit(noise).Rdata")

## Extracting AIC values from fitted models.
### tryCatch() function outputs 'NA' in case of any error.
All_AIC <- t(sapply(ParNames, function(x) sapply(Aval_Dist, function(z) tryCatch(AIC(fit[[x]][[z]]), error = function(e) NA))))
## Finding the PDF with the minimum AIC
Min_AIC <- apply(All_AIC, MARGIN = 1, function(x) names(which.min(x)))
## Saving the AIC and the best fitted model as a '.csv' file.
write.csv(cbind(All_AIC, Min_AIC), file = "AIC(noise).csv")

rm(i,j,ParNames,Aval_Dist,fit,All_AIC,Min_AIC)

#########################
##Proportion parameters##
#########################
# Finding parameters of proportion type.
ParNames <- rownames(D501)[D501$DataType == "Proportion"]
## Chosen probability distribution functions (PDFs) developed for discrete data ranging between [0,1).
Aval_Dist <- c("BB",	"BI",	"ZABB",	"ZABI","ZIBB","ZIBI")
### An object to save fitting results.
fit <- NULL

for (i in ParNames) {
  for (j in Aval_Dist) {
    ### A data frame of proportion values.
    #### y: A vector of (non-negative integer) quantiles.
    #### N: A vector of binomial denominators.
    temp <- NULL
    temp <- data.frame(y = Diploid_WT$n[,i], N = Diploid_WT$N[,i])
    ### Fitting various distributions developed for discrete data ranging between [0,1).
    ### tryCatch() function outputs 'NA' in case of any error.
    fit[[i]][[j]] <-  tryCatch(gamlss(cbind(y, N-y) ~ 1, data = na.omit(temp), family = j, control = gcont), error = function(e) NA)
  }
}
save(fit, file = "fit(proprtion).Rdata")

## Extracting AIC values from fitted models.
### tryCatch() function outputs 'NA' in case of any error.
All_AIC <- t(sapply(ParNames, function(x) sapply(Aval_Dist, function(z) tryCatch(AIC(fit[[x]][[z]]), error = function(e) NA))))
## Finding the PDF with the minimum AIC
Min_AIC <- apply(All_AIC, MARGIN = 1, function(x) names(which.min(x)))
## Saving the AIC and the best fitted model as a '.csv' file.
write.csv(cbind(All_AIC, Min_AIC), file = "AIC(proprtion).csv")

rm(i,j,ParNames,Aval_Dist,fit,All_AIC,Min_AIC)

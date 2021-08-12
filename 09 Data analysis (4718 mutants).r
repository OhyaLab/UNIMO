# To test the effectiveness of the proposed pipeline, we reanalyzed morphological variations of
# the 4718 yeast nonessential gene mutants by comparing them to a data set of haploid wild-type yeast
# strains (his3; 109 replicates) in the manner of the generalize linear model. 

# Models of the probability distributions for the 490 unimodal CalMorph parameters were determined 
# (see '04 Probability distribution function (114 WT).r' and '06 Modality  (114 WT).r' files)
# to accommodate the pre-defined statistical models.
# Then, we calculated the P value (two-sided test) as the deviation of each mutant from the wild-type
# (WT; null model) using functions in the 'gamlss' package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46).
# False discovery rate (FDR), the rate of type I errors in the rejected null hypothesis due to
# multiple comparisons, was estimated by the 'qvalue' package (Storey, 2002; J. R. Stat. Soc. B, 64(3): 479-498). 

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

## 'BiocManager' package is used to install 'qvalue' package.
### Looking for the 'BiocManager' package among the installed ones.
if("BiocManager" %in% rownames(installed.packages())){
  ### Loading 'BiocManager' package into R environment.
  library(BiocManager)
} else {install.packages("BiocManager")}

## 'qvalue' package is used to estimate FDR.
### Looking for the 'qvalue' package among the installed ones.
if("qvalue" %in% rownames(installed.packages())){
  ### Loading 'qvalue' package into R environment.
  library(qvalue)
  ### If 'qvalue' package is not installed, installing and loading it.
} else {BiocManager::install("qvalue"); library(qvalue)}


## Loading morphological data of 109 replicates of WT strains + 4718 yeast mutants.
### Note: first 109 rows are WT replicates.
### This data is already transformed and prepared in '08 Data preparation (4718 mutants).r' file.
load("NonEssential.Rdata")

# pBB() function of gamlss package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46)
# defines probability mass function of the beta-binomial distribution at any given sample
# between zero and one (0 < y < 1). In case of very large quantiles (numerator equals denominator)
# very small p values (less than 1E-13) for upper-tail test fail to be estimated.
# Following function (pBB2) solves this problem.
source("pBB2.R")

## Information of the CalMorph 501 morphological parameters.
## This information is needed to find the type of the data under 'DataType' column.
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)
## Finding unimodal parameters under 'Modality' column.
UniModal <- rownames(D501)[D501$Modality == 1]

## Controlling 'gamlss' fitting:
### Arguments are explained in gamlss.control() function of 'gamlss' package [help(gamlss.control, package = "gamlss")]:
#### n.cyc: The number of cycles of the algorithm.
#### trance: Logical; whether to print at each iteration.
gcont <- gamlss.control(n.cyc = 200L, trace = TRUE)

# Logistic function to accommodate 'logit' link function
logistic <- function(x) 1/(1+exp(-x))
#############################################
# Fitting a generalized linear model based on the predefined probability distributions 
# (see '04 Probability distribution function (114 WT).r' file).
# Here, we used 109 replicate of WT strain to defined the null model.

## An object to save the fitting results
glmres <-  list()

for (i in rownames(D501)) {
  ## Objects to defined the data and save the fitting result in each iteration
  temp <- fit <- NULL
  ###############
  # 1) Non-negative parameters (183 CalMorph parameters)
  ## 1-1) Gamma distribution (177 CalMorph parameters)
  if(D501[i,"PDF"] == "Gamma"){
    ### A data frame of:
    #### y: Non-negative values of 109 WT replicates as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ### Gamma distribution accepts data range of [0, +Inf).
    temp <- temp[0 < temp$y,,drop = FALSE]
    ## Fitting a generalized linear model give a Gamma distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = GA, control = gcont), error = function(e) NA)
  }
  ## 1-2) Weibull distribution (6 CalMorph parameters)
  if(D501[i,"PDF"] == "Weibull"){
    ### A data frame of:
    #### y: Non-negative values of 109 WT replicates as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ### Weibull distribution accepts data range of [0, +Inf).
    temp <- temp[0 < temp$y,,drop = FALSE]
    ## Fitting a generalized linear model give a Weibull distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = WEI, control = gcont), error = function(e) NA)
  }
  ###############
  ## 2) Ratio parameters (37 CalMorph parameters)
  ## 2-1) Beta distribution (36 CalMorph parameters)
  if(D501[i,"PDF"] == "Beta"){
    ### A data frame of:
    #### y: Ratio values of 109 WT replicates as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ### Beta distribution accepts data range of (0,1).
    temp <- temp[temp$y < 1,,drop = FALSE]
    ## Fitting a generalized linear model give a Beta distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = BE, control = gcont), error = function(e) NA)
  }
  ## 2-2) Logit normal distribution (1 CalMorph parameter)
  if(D501[i,"PDF"] == "Logit Normal"){
    ### A data frame of:
    #### y: Ratio values of 109 WT replicates as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ### Logit normal distribution accepts data range of (0,1).
    temp <- temp[temp$y < 1,,drop = FALSE]
    ## Fitting a generalized linear model give a Logit normal distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = LOGITNO, control = gcont), error = function(e) NA)
  }
  ###############
  # 3) Noise parameters (220 CalMorph parameters)
  ## 3-1) Gaussian distribution (209 CalMorph parameters)
  if(D501[i,"PDF"] == "Gaussian"){
    ### A data frame of:
    #### y: Noise values of 109 WT replicates as the dependent variable.
    ### Gaussian distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ## Fitting a generalized linear model give a Gaussian distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = NO, control = gcont), error = function(e) NA)
  }
  ## 3-2) Logistic distribution (5 CalMorph parameters)
  if(D501[i,"PDF"] == "Logistic"){
    ### A data frame of:
    #### y: Noise values of 109 WT replicates as the dependent variable.
    ### Logistic distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ## Fitting a generalized linear model give a Logistic distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = LO, control = gcont), error = function(e) NA)
  }
  ## 3-3) Reverse Gumbel distribution (6 CalMorph parameters)
  if(D501[i,"PDF"] == "Reverse Gumbel"){
    ### A data frame of:
    #### y: Noise values of 109 WT replicates as the dependent variable.
    ### Reverse Gumbel distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    ## Fitting a generalized linear model give a Reverse Gumbel distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = RG, control = gcont), error = function(e) NA)
  }
  ###############
  # 4) Proportion parameters (61 CalMorph parameters)
  ## 4-1) Binomial distribution (23 CalMorph parameters)
  if(D501[i,"PDF"] == "Binomial"){
    ### A data frame of:
    #### y: Vector of (non-negative integer) quantiles.
    #### N: Vector of binomial denominators.
    temp <- data.frame(y = AllNonEss$n[1:109,i], N = AllNonEss$N[1:109,i])
    ## Fitting a generalized linear model give a Binomial distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(cbind(y, N-y) ~ 1, data = na.omit(temp), family = BI, control = gcont), error = function(e) NA)
  }
  ## 4-2) Beta-Binomial distribution (38 CalMorph parameters) 
  if(D501[i,"PDF"] == "Beta-binomial"){ 
    ### A data frame of:
    #### y: Vector of (non-negative integer) quantiles.
    #### N: Vector of binomial denominators.
    temp <- data.frame(y = AllNonEss$n[1:109,i], N = AllNonEss$N[1:109,i])
    ## Fitting a generalized linear model give a Beta-Binomial distribution.
    ## tryCatch() function outputs 'NA' in case of any error.
    glmres[[i]] <- tryCatch(gamlss(cbind(y, N-y) ~ 1, data = na.omit(temp), family = BB, control = gcont), error = function(e) NA)
  }
}
rm(i)
###################################
# Estimating Q values given
## 1) predefined probability distribution of each parameter (see '04 Probability distribution function (114 WT).r' file)
## 2) unimodality of the parameter (see '06 Modality  (114 WT).r' files)
## 3) hyper-parameters of the null distribution (defined in lines 79-184)
## 4) two-sided test.

# An object to save P-values
PValue <- matrix(NA, nrow = nrow(AllNonEss$Trans), ncol = ncol(AllNonEss$Trans), dimnames = list(rownames(AllNonEss$Trans), colnames(AllNonEss$Trans)))


for (i in rownames(D501)) {
  ## Objects to defined the data and P values of the lower- and upper-tail in each iteration.
  temp <- res<- NULL
  ###############
  # 1) Non-negative parameters (183 CalMorph parameters)
  ## 1-1) Gamma distribution (177 CalMorph parameters)
  if(D501[i,"PDF"] == "Gamma"){
    ### A data frame of:
    #### y: Non-negative values as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### Gamma distribution accepts data range of [0, +Inf).
    temp <- temp[0 < temp$y,,drop = FALSE]
    ### EStimating P value given the lower-tail of a Gamma distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <- pGA(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Gamma distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <-  pGA(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
  }
  ## 1-2) Weibull distribution (6 CalMorph parameters)
  if(D501[i,"PDF"] == "Weibull"){
    ### A data frame of:
    #### y: Non-negative values as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### Weibull distribution accepts data range of [0, +Inf).
    temp <- temp[0 < temp$y,,drop = FALSE]
    ### EStimating P value given the lower-tail of a Weibull distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <- pWEI(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Weibull distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <- pWEI(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
  }
  ###############
  ## 2) Ratio parameters (37 CalMorph parameters)
  ## 2-1) Beta distribution (36 CalMorph parameters)
  if(D501[i,"PDF"] == "Beta"){
    ### A data frame of:
    #### y: Ratio values as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### Beta distribution accepts data range of (0,1).
    temp <- temp[temp$y < 1,,drop = FALSE]
    ### EStimating P value given the lower-tail of a Beta distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <- pBE(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = logistic(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Beta distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <- pBE(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = logistic(glmres[[i]]$sigma.coefficients), lower.tail = F)
  }
  ## 2-2) Logit normal distribution (1 CalMorph parameter)
  if(D501[i,"PDF"] == "Logit Normal"){
    ### A data frame of:
    #### y: Ratio values as the dependent variable.
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### Logit normal distribution accepts data range of (0,1).
    temp <- temp[temp$y < 1,,drop = FALSE]
    ### EStimating P value given the lower-tail of a Logit Normal distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <- pLOGITNO(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Logit Normal distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <- pLOGITNO(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
  }
  ###############
  # 3) Noise parameters (220 CalMorph parameters)
  ## 3-1) Gaussian distribution (209 CalMorph parameters)
  if(D501[i,"PDF"] == "Gaussian"){
    ### A data frame of:
    #### y: Noise values as the dependent variable.
    ### Gaussian distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### EStimating P value given the lower-tail of a Gaussian distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <- pNO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Gaussian distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <- pNO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
  }
  ## 3-2) Logistic distribution (5 CalMorph parameters)
  if(D501[i,"PDF"] == "Logistic"){
    ### A data frame of:
    #### y: Noise values as the dependent variable.
    ### Logistic distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### EStimating P value given the lower-tail of a Logistic distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <-  pLO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Logistic distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <- pLO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
  }
  ## 3-3) Reverse Gumbel distribution (6 CalMorph parameters)
  if(D501[i,"PDF"] == "Reverse Gumbel"){
    ### A data frame of:
    #### y: Noise values as the dependent variable.
    ### Reverse Gumbel distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = AllNonEss$Trans[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### EStimating P value given the lower-tail of a Reverse Gumbel distributed data (one-sided test).
    res[rownames(temp),"LowerTail"] <-  pRG(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    ### EStimating P value given the upper-tail of a Reverse Gumbel distributed data (one-sided test).
    res[rownames(temp),"UpperTail"] <- pRG(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    
  }
  ###############
  # 4) Proportion parameters (61 CalMorph parameters)
  ## 4-1) Binomial distribution (23 CalMorph parameters)
  if(D501[i,"PDF"] == "Binomial"){
    ### Weight for those replicates with zero or one values
    w <- ifelse(AllNonEss$n[,i] == 0 | AllNonEss$n[,i] == AllNonEss$N[,i], FALSE, TRUE)
    ### A data frame of:
    #### n: Vector of (non-negative integer) quantiles.
    #### N: Vector of binomial denominators.
    temp <- data.frame(n = AllNonEss$n[,i], N = AllNonEss$N[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### Observations should be passed to the pBI() function one by one.
    for (j in rownames(temp)[w]) {
      ### EStimating P value given the lower-tail of a Binomial distributed data (one-sided test).
      res[j,"LowerTail"] <- pBI(q = temp[j,"n"], bd = temp[j,"N"], mu = logistic(glmres[[i]]$mu.coefficients), lower.tail = T)
      ### EStimating P value given the upper-tail of a Binomial distributed data (one-sided test).
      res[j,"UpperTail"] <- pBI(q = temp[j,"n"], bd = temp[j,"N"], mu = logistic(glmres[[i]]$mu.coefficients), lower.tail = F)
    }
    res[is.na(res)] <- 1
  }
  ## 4-2) Beta-Binomial distribution (38 CalMorph parameters) 
  if(D501[i,"PDF"] == "Beta-binomial"){
    ### Weight for those replicates with zero or one values
    w <- ifelse(AllNonEss$n[,i] == 0 | AllNonEss$n[,i] == AllNonEss$N[,i], FALSE, TRUE)
    ### A data frame of:
    #### n: Vector of (non-negative integer) quantiles.
    #### N: Vector of binomial denominators.
    temp <- data.frame(n = AllNonEss$n[,i], N = AllNonEss$N[,i])
    ### A matrix to save P values of lower- and upper-tail
    res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(AllNonEss$Trans), c("LowerTail","UpperTail")))
    ### Observations should be passed to the pBB2() function one by one.
    for (j in rownames(temp)[w]) {
      ### EStimating P value given the lower-tail of a Beta-Binomial distributed data (one-sided test).
      res[j,"LowerTail"] <-  pBB2(q = temp[j,"n"], bd = temp[j,"N"], mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
      ### EStimating P value given the upper-tail of a Beta-Binomial distributed data (one-sided test).
      res[j,"UpperTail"] <-  pBB2(q = temp[j,"n"], bd = temp[j,"N"], mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    }
    res[is.na(res)] <- 1
  }
  ### Two-sided test: 2 * P value of a one-sided test.
  ### 2 * smaller P value between lower- and upper-tail.
  P2Sided <- 2 * (apply(res, MARGIN = 1, min))
  PValue[names(P2Side),i] <- P2Sided
  rm(P2Sided)
}
# Saving the estimated P values
save(PValue, file = "NonEssentialMut_PValue.Rdata")


# Making a a data frame of the P values
ThisStudyres <- do.call("cbind", PValue)
## qvalue() function of 'qvalue' package accepts data range of [0,1). 
ThisStudyres[ThisStudyres > 1] <- 1
## Finding unimodal parameters of 4718 mutants
Pvalue_res <- ThisStudyres[110:nrow(ThisStudyres),UniModal]
## Estimating Q values (FDR): No. of abnormal mutants
### FDR threshold
FDR <- .01
sum(rowSums(qvalue(Pvalue_res)$qvalues < FDR) > 0)

# We used mixture-model-based clustering, a flexible parametric method (Xu et al., 2014; Journal of Data Science 12(1): 175-196),
# to check the modality of the 501 CalMorph parameters (Ohya et al, 2005; PNAS 102(52): 19015-19020)
# according to their pre-defined distributions (defined in '04 Probability distribution function (114 WT).r' file).
# We used the 'gamlss.mx' package (Stasinopoulos et al., 2016), for all probability distributions except Gaussian.
# For the Gaussian-distributed parameters, we used the 'mclust' package given that univariate data have 
# equal volumes (Scrucca et al., 2016; R J. 8(1): 289).
# In all cases, we used the Bayesian information criterion (BIC) to compare the primary probability models
# that differed in the number of components:
#  BIC = ln(n) * k - 2ln(L)
## n: Number of observations (i.e., sample size).
## k: Number of hyper-parameters in the model multiplied by the number of the clusters (i.e., 1 <= c <= 10).
## L: Estimated maximum value of the likelihood function for the model.

# In the case of edge peak distributions, a model-free outlier-detection method (one-percentile deviation rule)
# was used to remove expected sets of outliers from the multimodal parameters.

# Authors: Ghanegolmohammadi et al. (2021)
##################################################
##################################################
##################################################
## Setting the working directory (not necessary to define)
setwd("...")
## Setting the initial seed (not necessary to set)
set.seed(123)

# Prerequisite
## 'gamlss.mx' package is used for mixture modeling
### Looking for the 'gamlss.mx' package among the installed ones
if("gamlss.mx" %in% rownames(installed.packages())){
  ### Loading 'gamlss.mx' package into R environment 
  library(gamlss.mx)
  ### If 'gamlss' package is not installed, installing and loading it.
} else {install.packages("gamlss.mx"); library(gamlss.mx)}

## 'mclust' package is used for Gaussian mixture modeling (GMM)
### Looking for the 'mclust' package among the installed ones
if("mclust" %in% rownames(installed.packages())){
  ### Loading 'mclust' package into R environment 
  library(mclust)
  ### If 'mclust' package is not installed, installing and loading it.
} else {install.packages("mclust"); library(mclust)}

## Loading 114 diploid Wild-type strains
### This data is already transformed and prepared in '02 Data preparation (114 WT).r' file.
load("Diploid_WT.Rdata")

## Information of the CalMorph 501 morphological parameters.
## This information is needed to find the type of the data under 'DataType' column and
## defining no. of hyper-parameters for each distribution under 'gamlss.function' column.
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

## Controlling 'gamlss' fitting:
### Arguments are explained in gamlss.control() function of 'gamlss' package [help(gamlss.control, package = "gamlss")]:
#### n.cyc: The number of cycles of the algorithm.
#### trance: Logical; whether to print at each iteration.
gcont <- gamlss.control(n.cyc = 200L, trace = TRUE)

## Controlling the gamlssMX function:
### Arguments are explained in MX.control() function of 'gamlss.mx' package [help(MX.control, package = "gamlss.mx")]:
#### cc: Convergent criterion for the EXpectation Maximization (EM).
#### n.cyc: The number of cycles for EM.
#### trance: Logical; whether to print at each iteration.
#### seed: A number for setting the seeds for starting values.
#### plot: Logical; whether to plot the sequence of global deviance up to convergence.
MXCon <- MX.control(cc = 1e-04, n.cyc = 200, trace = FALSE, 
                    seed = 123, plot = FALSE)
###################################
# Defining no. of hyper-parameters for calculating BIC
## Finding assigned distributions (defined in '04 Probability distribution function (114 WT).r' file).
Aval_Dist <- levels(D501[,"gamlss.function"])
## An object to save no. of hyper parameters in each distribution.
HypPar <- list()
## There are maximum four link functions in each distribution:
### mu: Compared to location hyper-parameter.
### sigma: Compared to dispersion hyper-parameter.
### nu: Compared to skewness hyper-parameter.
### tau: Compared to kurtosis hyper-parameter.
link_function <- c("mu","sigma","nu","tau")
for (i in Aval_Dist) {
  ## Getting head of the function coded in 'gamlss' package.
  temp <- NULL
  temp <- head(get(i))
  ### extracting link functions from the text
  HypPar[[as.character(i)]] <- gsub('[\\"\\"]', "", unlist(regmatches(temp, gregexpr('\\".*?\\"', temp))[1:2]))
  names(HypPar[[i]]) <- link_function[1: length(HypPar[[as.character(i)]])]
}
rm(i,temp)

## Defining max no. of clusters (k).
rep <- 10
#An object to save the results of mixture modeling of each parameters for K=1:rep.
MMAll <- list()
#An object to save BIC values.
BICAll <- list()

for(i in rownames(D501)){
  for(j in 1:rep){
    ## An object to save the mixture modeling for 1 <= c <= 10
    Modality.fit <- NULL
    ## An Object to defined the data
    temp <- NULL
    ###############
    # 1) Non-negative parameters (183 CalMorph parameters)
    ## 1-1) Gamma distribution (177 CalMorph parameters)
    if(D501[i,"PDF"] == "Gamma"){
      ### A data frame of:
      #### y: Non-negative values as the dependent variable.
      temp <- data.frame(y = His3Dip$Trans[,i])
      ### Gamma distribution accepts data range of [0, +Inf).
      temp <- temp[0 < temp$y,,drop = FALSE]
      ## Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Gamma distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(y ~ 1, data = na.omit(temp), family = GA, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ## 1-2) Weibull distribution (6 CalMorph parameters)
    if(D501[i,"PDF"] == "Weibull"){
      ### A data frame of:
      #### y: Non-negative values as the dependent variable.
      temp <- data.frame(y = His3Dip$Trans[,i])
      ### Weibull distribution accepts data range of [0, +Inf).
      temp <- temp[0 < temp$y,,drop = FALSE]
      ## Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Weibull distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(y ~ 1, data = na.omit(temp), family = WEI, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ###############
    ## 2) Ratio parameters (37 CalMorph parameters)
    ## 2-1) Beta distribution (36 CalMorph parameters)
    if(D501[i,"PDF"] == "Beta"){
      ### A data frame of:
      #### y: Ratio values as the dependent variable.
      temp <- data.frame(y = His3Dip$Trans[,i])
      ### Beta distribution accepts data range of (0,1).
      temp <- temp[temp$y < 1,,drop = FALSE]
      ## Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Beta distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(y ~ 1, data = na.omit(temp), family = BE, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ## 2-2) Logit normal distribution (1 CalMorph parameter)
    if(D501[i,"PDF"] == "Logit Normal"){
      ### A data frame of:
      #### y: Ratio values as the dependent variable.
      temp <- data.frame(y = His3Dip$Trans[,i])
      ### Logit normal distribution accepts data range of (0,1).
      temp <- temp[temp$y < 1,,drop = FALSE]
      ## Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Logit normal distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family = LOGITNO, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ###############
    # 3) Noise parameters (220 CalMorph parameters)
    ## 3-1) Gaussian distribution (209 CalMorph parameters)
    if(D501[i,"PDF"] == "Gaussian") {
      ### A data frame of:
      #### y: Noise values as the dependent variable.
      ### Gaussian distribution accepts data range of (-Inf,+Inf).
      temp <- data.frame(y = His3Dip$Trans[,i])
      ### Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Gaussian distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(Mclust(temp, G = j, modelNames = "E"), error = function(e) NA)
    }
    ## 3-2) Logistic distribution (5 CalMorph parameters)
    if(D501[i,"PDF"] == "Logistic") {
      ### A data frame of:
      #### y: Noise values as the dependent variable.
      ### Logistic distribution accepts data range of (-Inf,+Inf).
      temp <- data.frame(y = His3Dip$Trans[,i])
      ## Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Logistic distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(y ~ 1, data = na.omit(temp), family = LO, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ## 3-3) Reverse Gumbel distribution (6 CalMorph parameters)
    if(D501[i,"PDF"] == "Reverse Gumbel") {
      ### A data frame of:
      #### y: Noise values as the dependent variable.
      ### Reverse Gumbel distribution accepts data range of (-Inf,+Inf).
      temp <- data.frame(y = His3Dip$Trans[,i])
      ### Deleting one-percentile of the data as outliers.
      ### Run the following two lines for deleting the outliers.
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop = FALSE]
      ## Mixture model clustering, give a mixture of Reverse Gumbel distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(y ~ 1, data = na.omit(temp), family = RG, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ###############
    # 4) Proportion parameters (61 CalMorph parameters)
    ## 4-1) Binomial distribution (23 CalMorph parameters)
    if(D501[i,"PDF"] == "Binomial"){
      ### A data frame of:
      #### y: Proportion values as the dependent variable.
      temp.proportion <- data.frame(y = His3Dip$Trans[,i])
      ## Deleting one-percentile of the data as outliers.
      ### Run the following five lines for deleting the outliers.
      res <- quantile(temp.proportion$y, probs = seq(0,1,.005))
      NO_Outlier <- temp.proportion[temp.proportion$y > res["0.5%"] &  temp.proportion$y < res["99.5%"],,drop = FALSE]
      if(all(His3Dip$Trans[,i] == 0)){ # When all replicates are zero (e.g., "C123_C")
        ### A data frame of:
        #### n: Vector of (non-negative integer) quantiles.
        #### N: Vector of binomial denominators.
        temp <- data.frame(n = His3Dip$n[2:113,i], N = His3Dip$N[2:113,i])
      } else{
        ### A data frame of:
        #### n: Vector of (non-negative integer) quantiles.
        #### N: Vector of binomial denominators.
        temp <- data.frame(n = His3Dip$n[rownames(NO_Outlier),i], N = His3Dip$N[rownames(NO_Outlier),i])
      }
      ## Mixture model clustering, give a mixture of Binomial distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(cbind(n, N-n) ~ 1, data = na.omit(temp), family = BI, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    ## 4-2) Beta-Binomial distribution (38 CalMorph parameters) 
    if(D501[i,"PDF"] == "Beta-binomial") {
      ### A data frame of:
      #### y: Proportion values as the dependent variable.
      temp.proportion <- data.frame(y = His3Dip$Trans[,i])
      ## Deleting one-percentile of the data as outliers.
      ### Run the following five lines for deleting the outliers.
      res <- quantile(temp.proportion$y, probs = seq(0,1,.005))
      NO_Outlier <- temp.proportion[temp.proportion$y > res["0.5%"] &  temp.proportion$y < res["99.5%"],,drop = FALSE]
      if(all(His3Dip$Trans[,i] == 0)){ # When all replicates are zero
        ### A data frame of:
        #### n: Vector of (non-negative integer) quantiles.
        #### N: Vector of binomial denominators.
        temp <- data.frame(n = His3Dip$n[2:113,i], N = His3Dip$N[2:113,i])
      } else{
        ### A data frame of:
        #### n: Vector of (non-negative integer) quantiles.
        #### N: Vector of binomial denominators.
        temp <- data.frame(n = His3Dip$n[rownames(NO_Outlier),i], N = His3Dip$N[rownames(NO_Outlier),i])
      }
      ## Mixture model clustering, give a mixture of Beta-Binomial distributed data and 'j' clusters.
      ## tryCatch() function outputs 'NA' in case of any error.
      Modality.fit <- tryCatch(gamlssMX(cbind(n, N-n) ~ 1, data = na.omit(temp), family = BB, K = j,
                                        control = MXCon, g.control = gCont), error = function(e) NA)
    }
    MMAll[[i]][[paste("Cluster", j, sep = "=")]] <- Modality.fit
    
    # Calculating BIC: BIC=ln(n)k-2ln(L)
    ## k = The number of estimated hyper-parameters in the model * no. of clusters
    k <- length(HypPar[[D501[i,"gamlss.function"]]]) * j
    ## n = The number of observations, or equivalently (sample size)
    n <- nrow(temp)
    ## L = The maximum value of the likelihood function for the model
    if(D501[i,"PDF"] == "Gaussian"){
      L <- tryCatch(Modality.fit$loglik, error=function(e) NA)
    }else{
      L <- tryCatch(logLik(Modality.fit), error=function(e) NA)
    }
    ### BIC
    BICAll[[i]][[paste("Cluster", j, sep = "=")]] <- tryCatch((log(n)*k) - (2*as.numeric(L)), error=function(e) NA)
    rm(k,L,n)
  }
}
# Saving mixture models
save(MMAll, file = "MMAll.rdata")
# Saving BIC values
save(BICAll, file = "BICAll.rdata")
# Finding minimum BIC for 1<= k <= 10 of CalMorph parameters.
FinalBICALL <- do.call("rbind", BICAll)
Modality <- tryCatch(apply(FinalBICALL, MARGIN = 1, which.min), error=function(e) NA)
# Models of the probability distributions for each of the 501 CalMoprh parameters (Ohya et al, 2005; PNAS 102(52): 19015-19020)
# are needed to be determined to accommodate the statistical model used in the generalized linear modeling (GLM).
# CalMorph generates a variety of measurements that can be categorized into four groups. For each 
# group (i.e., data type), there is a well-known distribution that conventional has been used by statisticians.
# These common distributions are called reference distributions in this study. 

## 1) Non-negative continuous values [0, +Inf) includes 183 CalMorph parameters: Gamma (GA) distribution
## 2) Bounded continuous values; i.e., ratio [0,1] includes 37 CalMorph parameters: Beta (BE) distribution
## 3) Real continuous values; i.e., noise (-Inf, +Inf) includes 220 CalMorph parameters: Gaussian (NO) distribution
## 4) Finite discrete values; i.e., proportion [0,1) includes 61 CalMorph parameters: Binomial (BI) distribution

# Here, we check the equality of each parameter with the presumed reference distribution by considering
# significant deviation from related theoretical distribution using 
## 1) Kolmogorov-Smirnov test [help(ks.test, package = "stats")] for
### 1-1) Non-negative continuous values [0, +Inf).
### 1-2) Bounded continuous values; i.e., ratio [0,1].
### 1-3) Finite discrete values; i.e., proportion [0,1).
## 2) Shapiro-Wilk test of normality [help(shapiro.test, package = "stats")] for real continuous values; i.e., noise (-Inf, +Inf).

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
### This data is already transformed and prepared in '02 Data preparation (114 WT).r' file.
load("Diploid_WT.Rdata")

## Information of the CalMorph 501 morphological parameters.
## This information is needed to find the type of the data under 'DataType' column.
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

## Controlling 'gamlss' fitting:
### Arguments are explained in gamlss.control() function of 'gamlss' package [help(gamlss.control, package = "gamlss")]:
#### n.cyc: The number of cycles of the algorithm.
#### trance: Logical; whether to print at each iteration.
gcont <- gamlss.control(n.cyc = 200L, trace = TRUE)

# Logistic function to accommodate 'logit' link function
logistic <- function(x) 1/(1+exp(-x))
###################################
# Objects to save outputs of Kolmogorov-Smirnov and Shapiro-Wilk tests.
Kolmogorov <- list()
Shapiro <- list()

for (i in rownames(D501)) {
  ## Non-negative parameters.
  if(D501[i, DataType]== "Non-negative"){
    ### A data frame of non-negative values.
    temp <- NULL
    temp <- data.frame(y = Diploid_WT$Trans[,i])
    ### Fitting a Gamma distribution.
    ### tryCatch() function outputs 'NA' in case of any error.
    fit <- NULL
    fit <-  tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = GA, control = gcont), error = function(e) NA)
    ### Kolmogorov-Smirnov test given a theoretical Gamma distribution
    ### and estimated 'mu' and 'sigma' hyper-parameters.
    #### mu: Link function of the 'mu' hyper-parameter of Gamma distribution is 'log';
    #### thus, it should be exp() transformed; see help(GA, package = "gamlss.dist").
    #### sigma: Link function of the 'sigma' hyper-parameter of Gamma distribution is 'log';
    #### thus, it should be exp() transformed; see help(GA, package = "gamlss.dist").
    Kolmogorov[[i]] <- ks.test(Diploid_WT$Trans[,i], "pGA", mu = exp(fit$mu.coefficients),
                               sigma = exp(fit$sigma.coefficients), alternative = "two.sided")
  }
  ## Ratio parameters
  if(D501[i, DataType]== "Ratio"){
    ### A data frame of ratio values.
    temp <- NULL
    temp <- data.frame(y = Diploid_WT$Trans[,i])
    ### Fitting a Beta distribution.
    ### tryCatch() function outputs 'NA' in case of any error.
    fit <- NULL
    fit <- NULL
    fit <-  tryCatch(gamlss(y ~ 1, data = na.omit(temp), family = BE, control = gcont), error = function(e) NA)
    ### Kolmogorov-Smirnov test given a theoretical Beta distribution
    ### and estimated 'mu' and 'sigma' hyper-parameters.
    #### mu: Link function of the 'mu' hyper-parameter of Beta distribution is 'logit';
    #### thus, it should be logistic() transformed; see help(BE, package = "gamlss.dist").
    #### sigma: Link function of the 'sigma' hyper-parameter of Beta distribution is 'logit';
    #### thus, it should be logistic() transformed; see help(BE, package = "gamlss.dist").
    Kolmogorov[[i]] <- ks.test(Diploid_WT$Trans[,i], "pBE", mu = logistic(fit$mu.coefficients),
                               sigma = logistic(fit$sigma.coefficients), alternative = "two.sided")
  }
  ## Noise parameters
  if(D501[i, DataType]== "Noise"){
    ### Shapiro-Wilk test given a theoretical GAussian distribution.
    ### x: A numeric vector of a noise parameter.
    Shapiro[[i]] <- shapiro.test(x = Diploid_WT$Trans[,i])
  }
  ## Proportion parameters
  if(D501[i, DataType]== "Proportion"){
    ### A data frame of proportion values.
    #### y: A vector of (non-negative integer) quantiles.
    #### N: A vector of binomial denominators.
    temp <- NULL
    temp <- data.frame(y = Diploid_WT$n[,i], N = Diploid_WT$N[,i])
    ### Fitting a Binomial distribution.
    ### tryCatch() function outputs 'NA' in case of any error.
    fit <- NULL
    fit <-  tryCatch(gamlss(cbind(y, N-y) ~ 1, data = na.omit(temp), family = BI, control = gcont), error = function(e) NA)
    ### Kolmogorov-Smirnov test given a theoretical Binomial distribution
    ### randomly generated by 100-fold parametric bootstraps and estimated 'mu' hyper-parameter.
    #### mu: Link function of the 'mu' hyper-parameter of Binomial distribution is 'logit';
    #### thus, it should be logistic() transformed; see help(BI, package = "gamlss.dist").
    y <- table(His3Dip$N[,i])
    y <- data.frame(y)
    y$Var1 <- as.numeric(as.character(y$Var1))
    y <- data.frame(N = rep(y$Var1, y$Freq*100),
                    n = unlist(apply(as.matrix(y), 1, function(x) rBI(n = x[2]*100, bd = x[1], mu = logistic(fit$mu.coefficients)))))
    Kolmogorov[[i]] <- ks.test(x=finaldata$wt114$data[,i], y=y$n/y$N, alternative = "two.sided")
  }
}

# Saving the resutls of Kolmogorov-Smirnov and Shapiro-Wilk tests.
save(Kolmogorov, file = "Kolmogorov.Rdata")
save(Shapiro, file = "Shapiro.Rdata")

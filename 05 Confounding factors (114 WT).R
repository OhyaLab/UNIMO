# To minimize the effects of experimental error among replicates, a group of five confounding factors
# was considered, including a combination of different fluorescence filters for microscopy and
# the period of image acquisition (Ohnuki and Ohya, 2018; PLoS Biol., 16(5): e2005130).
# A generalized linear model was introduced by constructing a linear model (ANOVA) of the confounding factors:
# f(y) = (Beta_1 * MS1) + (Beta_2 * MS2) + (Beta_3 * MS2a) + (Beta_4 * MS2b) + (Beta_5 * MS3)
## Note: This 'Beta' is the Greek letter and it differs from the Beta distribution.
## y: The fitted value.
## f: An appropriate link function.
## Beta: Fixed effect of each confounding factor (0,1).
## MS1: Microscope #1 before replacement of the fluorescence filter.
## MS2: Microscope #2 before replacement of the fluorescence filter.
## MS2a: Microscope #2 after replacement of the fluorescence filter over time.
## MS2b: Microscope #2 after replacement of the fluorescence filter.
## MS3: Microscope #3.

# For every 501 CalMorph parameter (Ohya et al, 2005; PNAS 102(52): 19015-19020)
# the best-fitted model between the confounding factor model (CFM) and the null model was selected using
# AIC given the predefined probability distribution functions (PDFs).
## PDFs are defined in '04 Probability distribution function (114 WT).r' file.

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

## Loading morphological data of 114 diploid Wild-type strains
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

###################################
# A data frame to save the outputs of the fitting CFM and Null model
Confact <- data.frame(row.names = rownames(D501),
                      DataType = rep(NA, nrow(D501)), PDF = rep(NA, nrow(D501)),
                      CFM = rep(NA, nrow(D501)), CFM_df = rep(NA, nrow(D501)),
                      NLM = rep(NA, nrow(D501)), NLM_df = rep(NA, nrow(D501)),
                      LM = rep(NA, nrow(D501)))
###################################
# Generalized Linear Modeling (GLM)
for (i in rownames(D501)) {
  ## Finding data type for each parameter.
  Confact[i,"DataType"] <- as.character(D50[i,"DataType"])
  ## Finding PDF for each parameters (defined in '04 Probability distribution function (114 WT).r' file).
  Confact[i,"PDF"] <- as.character(D501[i,"PDF"])
  ## Objects to defined the data and save the fitting result
  res <- temp <- NULL
  
  ###############
  # 1) Non-negative parameters (183 CalMorph parameters)
  ## 1-1) Gamma distribution (177 CalMorph parameters)
  if(D501[i,"PDF"] == "Gamma"){
    ### A data frame of:
    #### y: Non-negative values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    temp <- data.frame(y = Diploid_WT$Trans[,i], x = Diploid_WT$info$mse)
    ### Gamma distribution accepts data range of [0, +Inf).
    temp <- temp[0 < temp$y,,drop=FALSE]
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = GA, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = GA, control=gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ## 1-2) Weibull distribution (6 CalMorph parameters)
  if(D501[i,"PDF"] == "Weibull"){
    ### A data frame of:
    #### y: Non-negative values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    temp <- data.frame(y = Diploid_WT$Trans[,i], x = Diploid_WT$info$mse)
    ### Weibull distribution accepts data range of [0, +Inf).
    temp <- temp[0 < temp$y,,drop=FALSE]
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = WEI, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = WEI, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ###############
  ## 2) Ratio parameters (37 CalMorph parameters)
  ## 2-1) Beta distribution (36 CalMorph parameters)
  if(D501[i,"PDF"] == "Beta"){
    ### A data frame of:
    #### y: Ratio values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    temp <- data.frame(y = Diploid_WT$Trans[,i], x = Diploid_WT$info$mse)
    ### Beta distribution accepts data range of (0,1).
    temp <- temp[temp$y < 1,,drop=FALSE]
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = BE, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = BE, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ## 2-2) Logit normal distribution (1 CalMorph parameter)
  if(D501[i,"PDF"] == "Logit Normal"){
    ### A data frame of:
    #### y: Ratio values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    temp <- data.frame(y = Diploid_WT$Trans[,i], x = Diploid_WT$info$mse)
    ### Logit normal distribution accepts data range of (0,1).
    temp <- temp[temp$y < 1,,drop=FALSE]
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = LOGITNO, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = LOGITNO, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ###############
  # 3) Noise parameters (220 CalMorph parameters)
  ## 3-1) Gaussian distribution (209 CalMorph parameters)
  if(D501[i,"PDF"] == "Gaussian") {
    ### A data frame of:
    #### y: Noise values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    ### Gaussian distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = Diploid_WT$Trans[,i], x = Diploid_WT$info$mse)
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = NO, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = NO, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ## 3-2) Logistic distribution (5 CalMorph parameters)
  if(D501[i,"PDF"] == "Logistic") {
    ### A data frame of:
    #### y: Noise values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    ### Logistic distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = LO, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = LO, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ## 3-3) Reverse Gumbel distribution (6 CalMorph parameters)
  if(D501[i,"PDF"] == "Reverse Gumbel") {
    ### A data frame of:
    #### y: Noise values as the dependent variable.
    #### x: Confounding factors as the explanatory variable.
    ### Reverse Gumbel distribution accepts data range of (-Inf,+Inf).
    temp <- data.frame(y = Diploid_WT$Trans[,i], x = Diploid_WT$info$mse)
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(y ~ x - 1, data = na.omit(temp), family = RG, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(y ~ 1, data = na.omit(temp), family = RG, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ###############
  # 4) Proportion parameters (61 CalMorph parameters)
  ## 4-1) Binomial distribution (23 CalMorph parameters)
  if(D501[i,"PDF"] == "Binomial"){
    ### A data frame of:
    #### n: Vector of (non-negative integer) quantiles.
    #### N: Vector of binomial denominators.
    #### x: Confounding factors as the explanatory variable.
    temp <- data.frame(n = His3Dip$n[,i], N = His3Dip$N[,i], x = Diploid_WT$info$mse)
    
    ### Model with confounding factors (CF).
    res$CFM <- gamlss(cbind(n, N-n) ~ x - 1, data = na.omit(temp), family = BI, control = gCont)
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- gamlss(cbind(n, N-n) ~ 1, data = na.omit(temp), family = BI, control = gCont)
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  ## 4-2) Beta-Binomial distribution (38 CalMorph parameters) 
  if(D501[i,"PDF"] == "Beta-binomial") {
    ### A data frame of:
    #### n: Vector of (non-negative integer) quantiles.
    #### N: Vector of binomial denominators.
    #### x: Confounding factors as the explanatory variable.
    temp <- data.frame(n = His3Dip$n[,i], N = His3Dip$N[,i], x = Diploid_WT$info$mse)
    
    ### Model with confounding factors (CF).
    res$CFM <- try(gamlss(cbind(n, N-n) ~ x - 1, data = na.omit(temp), family = BB, control = gCont))
    #### AIC of the fitted CF model.
    Confact[i,"CFM"] <- res$CFM$aic
    #### Degree of freedom of the fitted CF model.
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    
    ### Model without confounding factors (null model).
    res$NLM <- try(gamlss(cbind(n, N-n) ~ 1, data = na.omit(temp), family = BB, control = gCont))
    #### AIC of the fitted null model.
    Confact[i,"NLM"] <- res$NLM$aic
    #### Degree of freedom of the fitted null model.
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }
  
  ## Finding the best fitted linear model (LM) between CF and null models.
  Confact[i,"LM"] <- names(which.min(Confact[i,c("CFM","NLM")]))
}

# Saving the results.
write.csv(Confact, "ConfoundingFactorsModel.csv")
rm(i,res,fit,temp)
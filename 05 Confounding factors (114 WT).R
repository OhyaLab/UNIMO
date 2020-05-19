# Necessary codes to check fitness of null model vs confounding factor model
# for the 501 CalMorph parameters given the predefined PDFs

setwd("...")

set.seed(123)

# Prerequisites
if("gamlss" %in% rownames(installed.packages())){
  library(gamlss)
} else {install.packages("gamlss"); library(gamlss)}

## Loading 114 diploid Wild-type strains
load("Diploid_WT_Transformed.Rdata")

## Information of the CalMorph 501 morphological parameters
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

###################################
#Checking effects of confounding factors
Confact <- data.frame(row.names = rownames(D501),
                      DataType = rep(NA, nrow(D501)), PDF = rep(NA, nrow(D501)),
                      CFM = rep(NA, nrow(D501)), CFM_df = rep(NA, nrow(D501)),
                      NLM = rep(NA, nrow(D501)), NLM_df = rep(NA, nrow(D501)),
                      LM = rep(NA, nrow(D501)), Description = rep(NA, nrow(D501)))
###################################
#Controls for the gamlss fitting
gCont <- gamlss.control(n.cyc=200L, trace = FALSE)
#Generalized Linear Model (GLM)
for (i in rownames(D501)) {
  Confact[i,"DataType"] <- as.character(D501[i,"DataType"])
  Confact[i,"PDF"] <- as.character(D501[i,"PDF"])
  Confact[i,"Description"] <- as.character(D501[i,"Description"])
  res <- NULL
  temp <- NULL
  ###############
  if(D501[i,"PDF"] == "Gamma"){# Non-negative paramters (gamma distribution)
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    temp <- temp[0<temp$y,]
    # Model with confounding factors
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=GA, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=GA, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Weibull"){# Non-negative paramters (Weibull distribution)
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    temp <- temp[0<temp$y,]
    # Model with confounding factors
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=WEI, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=WEI, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Beta"){# Ratio paramters (beta distribution)
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    temp <- temp[temp$y<1,]
    # Model with confounding factors
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=BE, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=BE, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Logit Normal"){# Ratio paramters (logit Normal distribution)
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    temp <- temp[temp$y<1,]
    # Model with confounding factors
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=LOGITNO, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=LOGITNO, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Gaussian") {# Noise parameters (Gaussian distibution)
    # Model with confounding factors
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=NO, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=NO, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Logistic") {# Noise parameters (logistic distibution)
    # Model with confounding factors
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=LO, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=LO, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Reverse Gumbel") {# Noise parameters (reverse Gumbel distibution)
    # Model with confounding factors
    temp <- data.frame(y=Diploid_WT$Trans[,i], x=Diploid_WT$info$mse)
    res$CFM <- gamlss(y~x-1, data=na.omit(temp), family=RG, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(y~1, data=na.omit(temp), family=RG, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Binomial"){# Proportion parameters (binomial paramters)
    # Model with confounding factors
    temp <- data.frame(n=His3Dip$n[,i], N=His3Dip$N[,i], x=Diploid_WT$info$mse)
    res$CFM <- gamlss(cbind(n, N-n)~x-1, data=na.omit(temp), family=BI, control=gCont)
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- gamlss(cbind(n, N-n)~1, data=na.omit(temp), family=BI, control=gCont)
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
  }else if(D501[i,"PDF"] == "Beta-binomial") {# Proportion parameters (beta-binomial distibution)
    temp <- data.frame(n=His3Dip$n[,i], N=His3Dip$N[,i], x=Diploid_WT$info$mse)
    # Model with confounding factors
    res$CFM <- try(gamlss(cbind(n, N-n)~x-1, data=na.omit(temp), family=BB, control=gCont))
    Confact[i,"CFM"] <- res$CFM$aic
    Confact[i,"CFM_df"] <- res$CFM$df.fit
    # Model without confounding factors
    res$NLM <- try(gamlss(cbind(n, N-n)~1, data=na.omit(temp), family=BB, control=gCont))
    Confact[i,"NLM"] <- res$NLM$aic
    Confact[i,"NLM_df"] <- res$NLM$df.fit
    }
  Confact[i,"LM"] <- names(which.min(Confact[i,c("CFM","NLM")]))
}

write.csv(Confact, "ConfoundingFactorsModel.csv")
rm(i, fit, temp)
# Necessary codes to to check the modality of each 501 CalMorph parameters given the predefined PDFs

setwd("...")

set.seed(123)

# Prerequisites
if("gamlss.mx" %in% rownames(installed.packages())){
  library(gamlss.mx)
} else {install.packages("gamlss.mx"); library(gamlss.mx)}

if("mclust" %in% rownames(installed.packages())){
  library(mclust)
} else {install.packages("mclust"); library(mclust)}

## Loading 114 diploid Wild-type strains
load("Diploid_WT.Rdata")

## Information of the CalMorph 501 morphological parameters
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)
###################################################
###################################################
# Modality check
## Defining no. of hyper-parameters for calculating BIC
Aval_Dist <- levels(D501[,"gamlss.function"])
HypPar <- list()
link_function <- c("mu","sigma","nu","tau")
for (i in Aval_Dist) {
  temp <- NULL
  temp <- head(get(i))
  HypPar[[as.character(i)]] <- gsub('[\\"\\"]', "", unlist(regmatches(temp, gregexpr('\\".*?\\"', temp))[1:2]))
  names(HypPar[[i]]) <- link_function[1: length(HypPar[[as.character(i)]])]
}
rm(i,temp)

## Defining max no. of clusters (k)
rep <- 10
## Controls for the gamlssMX function
MXCon <- MX.control(cc = 1e-04, n.cyc = 200, trace = FALSE, 
                    seed = 123, plot = FALSE)
## Controls for the gamlss fitting
gCont <- gamlss.control(n.cyc=2000L, trace = FALSE)

#An onject to save the results of mixture modeling of each paramters for K=1:rep
MMAll <- list()
#An object to save BIC values
BICAll <- list()
#An object to save AIC values
AICAll <- list()

for(i in rownames(D501)){
  for(j in 1:rep){
    Modality.fit <- NULL
    temp <- NULL
    if(D501[i,"PDF"] == "Gamma"){# Non-negative paramters (Gamma distribution)
      temp <- data.frame(y=His3Dip$Trans[,i])
      temp <- temp[0<temp$y,,drop=FALSE]
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family=GA, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Weibull"){# Non-negative paramters (Weibull distribution)
      temp <- data.frame(y=His3Dip$Trans[,i])
      temp <- temp[0<temp$y,,drop=FALSE]
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family=WEI, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Beta"){# Ratio paramters (Beta distribution)
      temp <- data.frame(y=His3Dip$Trans[,i])
      temp <- temp[temp$y<1,,drop=FALSE]
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family=BE, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Logit Normal"){# Ratio paramters (Logit Normal distribution)
      temp <- data.frame(y=His3Dip$Trans[,i])
      temp <- temp[temp$y<1,,drop=FALSE]
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family=LOGITNO, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Gaussian") {# Gaussian distibution
      temp <- data.frame(y=His3Dip$Trans[,i])
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(Mclust(temp, G=j, modelNames = "E"), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Logistic") {# Logistic distibution
      temp <- data.frame(y=His3Dip$Trans[,i])
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family=LO, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Reverse Gumbel") {# Reverse Gumbel distibution
      temp <- data.frame(y=His3Dip$Trans[,i])
      # Run the following two lines for deleting the outliers
      res <- quantile(temp$y, probs = seq(0,1,.005))
      temp <- temp[temp$y > res["0.5%"] &  temp$y < res["99.5%"],,drop=FALSE]
      
      Modality.fit <- tryCatch(gamlssMX(y~1, data = na.omit(temp), family=RG, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Binomial"){# Binomial paramters
      # n: No. of success; N: Total no. of observations
      temp.proportion <- data.frame(y=His3Dip$Trans[,i])
      # Run the following 7 lines for deleting the outliers
      res <- quantile(temp.proportion$y, probs = seq(0,1,.005))
      NO_Outlier <- temp.proportion[temp.proportion$y > res["0.5%"] &  temp.proportion$y < res["99.5%"],,drop=FALSE]
      if(all(His3Dip$Trans[,i] == 0)){ # When all replicates are zero (e.g., "C123_C")
        temp <- data.frame(n=His3Dip$n[2:113,i], N=His3Dip$N[2:113,i])
      } else{
        temp <- data.frame(n=His3Dip$n[rownames(NO_Outlier),i], N=His3Dip$N[rownames(NO_Outlier),i])
      }
      
      Modality.fit <- tryCatch(gamlssMX(cbind(n, N-n)~1, data = na.omit(temp), family=BI, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }else if(D501[i,"PDF"] == "Beta-binomial") {# Beta-binomial distibution
      #n: No. of success; N: Total no. of observations
      temp.proportion <- data.frame(y=His3Dip$Trans[,i])
      # Run the following 7 lines for deleting the outliers
      res <- quantile(temp.proportion$y, probs = seq(0,1,.005))
      NO_Outlier <- temp.proportion[temp.proportion$y > res["0.5%"] &  temp.proportion$y < res["99.5%"],,drop=FALSE]
      if(all(His3Dip$Trans[,i] == 0)){ # When all replicates are zero
        temp <- data.frame(n=His3Dip$n[2:113,i], N=His3Dip$N[2:113,i])
      } else{
        temp <- data.frame(n=His3Dip$n[rownames(NO_Outlier),i], N=His3Dip$N[rownames(NO_Outlier),i])
      }
      
      Modality.fit <- tryCatch(gamlssMX(cbind(n, N-n)~1, data = na.omit(temp), family=BB, K=j,
                                        control = MXCon, g.control = gCont), error=function(e) NA)
    }
    MMAll[[i]][[paste("Cluster", j, sep = "=")]] <- Modality.fit
    
    # Calculating BIC: BIC=ln(n)k-2ln(L)
    # k = the number of estimated hyper-parameters in the model * no. of clusters
    k <- length(HypPar[[D501[i,"gamlss.function"]]]) * j
     #n is the number of observations, or equivalently (sample size)
    n <- nrow(temp)
    # L is the maximum value of the likelihood function for the model
    if(D501[i,"PDF"] == "Gaussian"){
      L <- tryCatch(Modality.fit$loglik, error=function(e) NA)
    }else{
      L <- tryCatch(logLik(Modality.fit), error=function(e) NA)
    }
    # BIC
    BICAll[[i]][[paste("Cluster", j, sep = "=")]] <- tryCatch((log(n)*k) - (2*as.numeric(L)), error=function(e) NA)
    rm(k,L,n)
  }
}

save(MMAll, file = "MMAll.rdata")
save(BICAll, file = "BICAll.rdata")

FinalBICALL <- do.call("rbind", BICAll)
Modality <- tryCatch(apply(FinalBICALL, MARGIN = 1, which.min), error=function(e) NA)
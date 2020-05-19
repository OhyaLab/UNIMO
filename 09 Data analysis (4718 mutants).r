# Necessary codes for finding significant morphological change of 4718 nonessential mutant give the new pipeline
setwd("...")

set.seed(123)

# Prerequisite
if("gamlss" %in% rownames(installed.packages())){
  library(gamlss)
} else {install.packages("gamlss"); library(gamlss)}

if("qvalue" %in% rownames(installed.packages())){
  library(qvalue)
} else {install.packages("qvalue"); library(qvalue)}

## Loading the transformed data
load("NonEssential.Rdata")

## Solving the problem for very small p values for upper-tail test of beta-binomial
source("pBB2.R")

## Supplementary table 4
D501 <- read.csv("501-InfoNew.csv", header = TRUE, row.name = 1)
UniModal <- rownames(D501)[D501$Modality == 1]
#############################################
# Fitting a distribution based on the previous results
## Controls for the gamlss fitting
gcont <- gamlss.control(n.cyc=200L, trace = FALSE)
## An onject to save the results
glmres <-  list()
for (i in rownames(D501)) {
  temp <- NULL
  fit <- NULL
  if(D501[i,"PDF"] == "Gamma"){ # Non-negative prameters (gamma distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    temp <- temp[0<temp$y,,drop=FALSE]
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=GA, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Weibull"){  # Non-negative prameters (Weibull distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    temp <- temp[0<temp$y,,drop=FALSE]
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=WEI, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Beta"){  # Ratio parameters (beta distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    temp <- temp[temp$y<1,,drop=FALSE]
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=BE, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Logit Normal"){  # Ratio parameters (logit normal distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    temp <- temp[temp$y<1,,drop=FALSE]
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=LOGITNO, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Gaussian"){  # Noise parameters (Gaussian distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=NO, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Logistic"){  # Noise parameters (logistic distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=LO, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Reverse Gumbel"){  # Noise parameters (reverse Gumbel distribution)
    temp <- data.frame(y = AllNonEss$Trans[1:109,i])
    glmres[[i]] <- tryCatch(gamlss(y~1, data=na.omit(temp), family=RG, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Binomial"){  # Proportion parameters (binomial distribution)
    temp <- data.frame(y = AllNonEss$n[1:109,i], N = AllNonEss$N[1:109,i])
    glmres[[i]] <- tryCatch(gamlss(cbind(y, N-y)~1, data=na.omit(temp), family=BI, control=gcont), error=function(e) NA)
  }
  if(D501[i,"PDF"] == "Beta-binomial"){  # Proportion parameters (beta-binomial distribution)
    temp <- data.frame(y = AllNonEss$n[1:109,i], N = AllNonEss$N[1:109,i])
    glmres[[i]] <- tryCatch(gamlss(cbind(y, N-y)~1, data=na.omit(temp), family=BB, control=gcont), error=function(e) NA)
  }
}
rm(i)
###################################
logistic <- function(x) 1/(1+exp(-x))
# An object to save P-values
PValue <- list()
# Calcualting P-values
for (i in rownames(D501)) {
  temp <- NULL
  res<- NULL
  if(D501[i,"PDF"] == "Gamma"){# Non-negative prameters (gamma distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    temp <- temp[0 < temp$y,,drop=FALSE]
    res$Low <- pGA(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <-  pGA(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Weibull"){# Non-negative prameters (Weibull distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    temp <- temp[0 < temp$y,,drop=FALSE]
    res$Low <- pWEI(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <- pWEI(q = temp$y, mu = exp(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Beta"){# Ratio prameters (beta distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    temp <- temp[temp$y < 1,,drop=FALSE]
    res$Low <- pBE(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = logistic(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <- pBE(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = logistic(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Logit Normal"){# Ratio prameters (logit Normal distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    temp <- temp[temp$y < 1,,drop=FALSE]
    res$Low <- pLOGITNO(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <- pLOGITNO(q = temp$y, mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Gaussian"){# Noise prameters (Gaussian distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    res$Low <- pNO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <- pNO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Logistic"){# Noise prameters (logistic distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    res$Low <-  pLO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <- pLO(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Reverse Gumbel"){# Noise prameters (reverse Gumbel distribution)
    temp <- data.frame(y = AllNonEss$Trans[,i])
    res$Low <-  pRG(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
    res$Up <- pRG(q = temp$y, mu = glmres[[i]]$mu.coefficients, sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Binomial"){# Noise prameters (binomial distribution)
    temp <- data.frame(n = AllNonEss$n[,i], N = AllNonEss$N[,i])
    for (j in 1:nrow(temp)) {
      res$Low[[j]] <- pBI(q = temp$n[j], bd = temp$N[j], mu = logistic(glmres[[i]]$mu.coefficients), lower.tail = T)
      res$Up[[j]] <- pBI(q = temp$n[j], bd = temp$N[j], mu = logistic(glmres[[i]]$mu.coefficients), lower.tail = F)
    }
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
  if(D501[i,"PDF"] == "Beta-binomial"){# Noise prameters (beta-binomial distribution)
    temp <- data.frame(n = AllNonEss$n[,i], N = AllNonEss$N[,i])
    for (j in 1:nrow(temp)) {
      res$Low[[j]] <-  pBB2(q = temp$n[j], bd = temp$N[j], mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = T)
      res$Up[[j]] <-  pBB2(q = temp$n[j], bd = temp$N[j], mu = logistic(glmres[[i]]$mu.coefficients), sigma = exp(glmres[[i]]$sigma.coefficients), lower.tail = F)
    }
    PValue[[i]] <- 2 * (apply(cbind(res$Low, res$Up), MARGIN = 1, min))
  }
}
save(PValue, file = "NonEssentialMut_PValue.Rdata")


# All P-values
ThisStudyres <- do.call("cbind", PValue)
## Dealing with the data out of range data
ThisStudyres[ThisStudyres > 1] <- 1
## Finding multimodal parameters and just mutantds
Pvalue_res <- ThisStudyres[110:nrow(ThisStudyres),UniModal]
## Estimateing Q values: No. of abnormal mutants
m <- .01
sum(rowSums(qvalue(Pvalue_res)$qvalues < m) > 0)

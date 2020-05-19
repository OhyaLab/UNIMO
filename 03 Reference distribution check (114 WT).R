# Necessary codes to check the equality of the data with the presumed reference distributions
# for the 501 CalMorph parameters given the data type

setwd("...")

set.seed(123)

# Prerequisites
if("gamlss" %in% rownames(installed.packages())){
  library(gamlss)
} else {install.packages("gamlss"); library(gamlss)}


## Loading 114 diploid Wild-type strains
load("Diploid_WT.Rdata")

## Information of the CalMorph 501 morphological parameters
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

gcont <- gamlss.control(n.cyc=200L, trace = TRUE)
logistic <- function(x) 1/(1+exp(-x))

###################################
# Finding each data type

ParNames <- rownames(D501)[D501$DataType == "Non-negative"]

Kolmogorov <- list()
Shapiro <- list()

for (i in rownames(D501)) {
  ## Non-negative parameters
  if(D501[i, DataType]== "Non-negative"){
    temp <- NULL
    temp <- data.frame(y=Diploid_WT$Trans[,i])
    fit <- NULL
    fit <-  tryCatch(gamlss(y~1, data=na.omit(temp), family=GA, control=gcont), error=function(e) NA)
    Kolmogorov[[i]] <- ks.test(Diploid_WT$Trans[,i], "pGA",
                               mu=exp(fit$mu.coefficients),
                               sigma=exp(fit$sigma.coefficients),
                               alternative = "two.sided")
  }
  ## Ratio parameters
  if(D501[i, DataType]== "Ratio"){
    temp <- NULL
    temp <- data.frame(y=Diploid_WT$Trans[,i])
    fit <- NULL
    fit <-  tryCatch(gamlss(y~1, data=na.omit(temp), family=BE, control=gcont), error=function(e) NA)
    Kolmogorov[[i]] <- ks.test(Diploid_WT$Trans[,i], "pBE",
                               mu=logistic(fit$mu.coefficients),
                               sigma=logistic(fit$sigma.coefficients),
                               alternative = "two.sided")
  }
  ## Noise parameters
  if(D501[i, DataType]== "Noise"){
    Shapiro[[i]] <- shapiro.test(Diploid_WT$Trans[,i])
  }
  ## Proportion parameters
  if(D501[i, DataType]== "Proportion"){
    temp <- NULL
    temp <- data.frame(y=Diploid_WT$n[,i], N=Diploid_WT$N[,i])
    fit <- NULL
    fit <-  tryCatch(gamlss(cbind(y, N-y)~1, data=na.omit(temp), family=BI, control=gcont), error=function(e) NA)
    y <- table(His3Dip$N[,i])
    y <- data.frame(y)
    y$Var1 <- as.numeric(as.character(y$Var1))
    y <- data.frame(N=rep(y$Var1, y$Freq*100),
                    n=unlist(apply(as.matrix(y), 1, function(x) rBI(n=x[2]*100, bd=x[1], mu=logistic(fit$mu.coefficients)))))
    Kolmogorov[[i]] <- ks.test(x=finaldata$wt114$data[,i], y=y$n/y$N, alternative = "two.sided")
  }
}

save(Kolmogorov, file = "Kolmogorov.Rdata")
save(Shapiro, file = "Shapiro.Rdata")

# Necessary codes to fit various probability distribution functions (PDFs)
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
###################################################
###################################################
# Finding parameters of each type
## Run the next two lines if you are analyzing non-negative parameters
ParNames <- rownames(D501)[D501$DataType == "Non-negative"]
## Chosen PDFs
Aval_Dist <- c("EXP","GA","GP","IG","IGAMMA","LOGNO","LOGNO2","PARETO1", "PARETO2","PARETO2o",
               "WEI","WEI2","WEI3","ZAGA","ZAIG")

## Run the next two lines if you are analyzing ratio parameters
ParNames <- rownames(D501)[D501$DataType == "Ratio"]
## Chosen PDFs
Aval_Dist <- c("BE","BEINF0","BEINF1","BEINF","BEOI","BEZI","LOGITNO","ZAGA")

## Run the next two lines if you are analyzing noise parameters
ParNames <- rownames(D501)[D501$DataType == "Noise"]
## Chosen PDFs
Aval_Dist <- c("GU","LO","NO","NO2","RG")

## Run the next two lines if you are analyzing proportion parameters
ParNames <- rownames(D501)[D501$DataType == "Proportion"]
## Chosen PDFs
Aval_Dist <- c("BB",	"BI",	"ZABB",	"ZABI","ZIBB","ZIBI")

###################################################
###################################################
# Fitting the chosen PDFs
## Controls for the gamlss fitting
gcont <- gamlss.control(n.cyc=200L, trace = TRUE)
## An object to save fitting results
fit <- list()
## Run following lines for non-negative OR ratio OR noise parameters
for (i in ParNames) {
  for (j in Aval_Dist) {
    temp <- NULL
    temp <- data.frame(y = Diploid_WT$Trans[,i])
    # Fitting distribution
    fit[[i]][[j]] <-  tryCatch(gamlss(y~1, data=na.omit(temp), family=j, control=gcont), error=function(e) list(NA))
  }
}
# ## Run following lines for proportion parameters
fit <- NULL
for (i in ParNames) {
  for (j in Aval_Dist) {
    temp <- NULL
    temp <- data.frame(y=Diploid_WT$n[,i], N=Diploid_WT$N[,i])
    #Fitting distribution
    fit[[i]][[j]] <-  tryCatch(gamlss(cbind(y, N-y)~1, data=na.omit(temp), family=j, control=gcont), error=function(e) NA)
  }
}
rm(i,j, temp)

save(fit, file = "BestDist.Rdata")

## Extracting AICs
All_AIC <- t(sapply(ParNames, function(x) sapply(Aval_Dist, function(z) tryCatch(AIC(fit[[x]][[z]]), error=function(e) NA))))

## Finding the PDF with the minimum AIC
Min_AIC <- apply(All_AIC, MARGIN = 1, function(x) names(which.min(x)))
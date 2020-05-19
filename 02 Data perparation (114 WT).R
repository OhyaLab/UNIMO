# Necessary codes for transformation of the 501 CalMorph parameters
# including:
# 1) Uncoupling the dependency between the coefficient of variation (CV) and the mean values (Levy $ Siegal, 2008; Yvert et al. 2013)
# 2) Handling inversion of some parameters

setwd("...")

set.seed(123)

# Prerequisite
## LOWESS regression function
source("lowess.fit.r")

## Information of the CalMorph 501 morphological parameters
D501 <- read.csv("CalMorph501.csv", header = TRUE, row.name = 1)

Diploid_WT <- NULL
## Loading 114 diploid Wild-type strains:
### Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114data.tsv
Diploid_WT$Data <- read.table("wt114data.tsv", header = TRUE, row.names = 1)
colnames(Diploid_WT$Data) <- chartr(".", "-", colnames(Diploid_WT$Data))
### Number of cells for proportion parameters: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114nmrt.tsv
Diploid_WT$n <- read.table("wt114nmrt.tsv", header = TRUE, row.names = 1)
### Number of cells in specimen for proportion parameters: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114dmnt.tsv
Diploid_WT$N <- read.table("wt114dmnt.tsv", header = TRUE, row.names = 1)
#############################################
# Data preparation
## An object to save the transformed values
Diploid_WT$Trans <- Diploid_WT$Data

## Decoupling of non-linearity for the noise parameters
idx <- rownames(D501)[grep("CV", rownames(D501))]
names(idx) <- idx
idx <- sub("CV", "", idx)

## Reading smooth span (defined by the lowest AIC)
### See "01 NoiseParameters.R"
load("SmoothSpan_WT114.rdata")
## Transformtation
for(i in names(idx)) {
  # x = Mean parameters
  # y = CV parameters
  temp <- lowess.fit(x = Diploid_WT$Data[,idx[i]], y = Diploid_WT$Data[,i],
                     ret = "all", f = SmoothSpan[i,1])
  # Noise = observed CV - predicted mean
  # Replacing the values of CV parameters by noise values
  Diploid_WT$Trans[,i] <- temp$res
}
rm(i,temp)

## Handling invers ("1/logit") and inver-double ("2/logit")
for(i in rownames(D501)) { 
  if(D501[i,"Link_function"] == "1/logit") {
    # Handling inverse
    Diploid_WT$Trans[,i] <- 1/Diploid_WT$Trans[,i]
  } else if(D501[i,"Link_function"] == "2/logit") {
    # Handling inverse-double (double reciprocal)
    Diploid_WT$Trans[,i] <- 2/Diploid_WT$Trans[,i]
  }
}
rm(i)

# Confounding factor information
load("Confounding_factors.Rdata")
Diploid_WT$info <- Confounding_factors

save(Diploid_WT, file = "Diploid_WT.Rdata")

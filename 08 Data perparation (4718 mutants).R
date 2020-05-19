# Necessary codes for transformation of the 501 CalMorph parameters
# including:
# 1) Uncoupling the dependency between the coefficient of variation (CV) and the mean values (Levy $ Siegal, 2008; Yvert et al. 2013)
# 2) Handling inversion of some parameters

setwd("...")

set.seed(123)

# Prerequisite
## LOWESS regression function
source("lowess.fit.r")

## Supplementary table 4
D501 <- read.csv("501-InfoNew.csv", header = TRUE, row.name = 1)

AllNonEss <- NULL
# Loading 4718 non-essential mutants:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718data.tsv
AllNonEss$Mut <- read.table("mt4718data.tsv", header = TRUE, row.names = 1)
colnames(AllNonEss$Mut) <- chartr(".","-", colnames(AllNonEss$Mut))
## Number of cells for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718nmrt.tsv
AllNonEss$Mut$n <- read.table("mt4718nmrt.tsv", header = TRUE, row.names = 1)
## Number of cells in specimen for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718dmnt.tsv
AllNonEss$Mut$N <- read.table("mt4718dmnt.tsv", header = TRUE, row.names = 1)

# Loading 122 haploid Wild-type strains:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122data.tsv
AllNonEss$WT <- read.table("wt122data.tsv", header = TRUE, row.names = 1)
colnames(AllNonEss$WT) <- chartr(".","-", colnames(AllNonEss$WT))
## Number of cells for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122nmrt.tsv
AllNonEss$WT$n <- read.table("wt122nmrt.tsv", header = TRUE, row.names = 1)
## Number of cells in specimen for proportion parameter: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122dmnt.tsv
AllNonEss$WT$N <- read.table("wt122dmnt.tsv", header = TRUE, row.names = 1)

### 13 replicates of the 122 WT strains were not used (see Suzuki et al., 2018, BMC Genomics)
res <- c("04his3-12", "04his3-13", "04his3-18", "04his3-22", "04his3-33", "04his3-54", 
         "04his3-55", "04his3-60", "04his3-70", "04his385", "his3cnt-1", "his3old-1", "his3old-5")
AllNonEss$WT <- AllNonEss$WT[!rownames(AllNonEss$WT) %in% res,]
AllNonEss$WT$n <- AllNonEss$WT$n[!rownames(AllNonEss$WT$n) %in% res,]
AllNonEss$WT$N <- AllNonEss$WT$N[!rownames(AllNonEss$WT$N) %in% res,]
rm(res)

### Final data
AllNonEss$Data <- rbind(AllNonEss$WT, AllNonEss$Mut)
AllNonEss$n <- rbind(AllNonEss$WT$n, AllNonEss$Mut$n)
AllNonEss$N <- rbind(AllNonEss$WT$N, AllNonEss$Mut$N)
#############################################
# Data preparation
## An object to save the transformed values
AllNonEss$Trans <- AllNonEss$Data

## Decoupling of non-linearity for the noise parameters
idx <- rownames(D501)[grep("CV", rownames(D501))]
names(idx) <- idx
idx <- sub("CV", "", idx)

## Reading smooth span (defined by the lowest AIC)
### See "01 NoiseParameters.R"
load("SmoothSpan_NonEssential.rdata")
## Transformtation
for(i in names(idx)) {
  # x = Mean parameters
  # y = CV parameters
  temp <- lowess.fit(x = AllNonEss$Data[,idx[i]], y = AllNonEss$Data[,i],
                     ret = "all", f = SmoothSpan[i,1])
  # Noise = observed CV - predicted mean
  # Replacing the values of CV parameters by noise values
  AllNonEss$Trans[,i] <- temp$res
}
rm(i,temp)

## Handling invers ("1/logit") and inver-double ("2/logit")
for(i in rownames(D501)) {
  if(D501[i,"Link_function"] == "1/logit") {
    # Handling inverse
    AllNonEss$Trans[,i] <- 1/AllNonEss$Trans[,i]
  } else if(D501[i,"Link_function"] == "2/logit") {
    # Handling inverse-double (double reciprocal)
    AllNonEss$Trans[,i] <- 2/AllNonEss$Trans[,i]
  }
}
rm(i)

save(AllNonEss, file = "NonEssential.Rdata")

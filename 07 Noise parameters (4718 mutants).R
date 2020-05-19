# Necessary code to calcualte the span for LOWESS regression

setwd("...")

set.seed(123)

# Prerequisite
if("gamlss" %in% rownames(installed.packages())){
  library(gamlss)
} else {install.packages("gamlss"); library(gamlss)}


AllHap <- NULL
# Loading 4718 non-essential mutants:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=mt4718data.tsv
AllHap$Mut <- read.table("mt4718data.tsv", header = TRUE, row.names = 1)
colnames(AllNonEss$Mut) <- chartr(".","-", colnames(AllNonEss$Mut))
# Loading 122 haploid Wild-type strains:
## Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt122data.tsv
AllHap$WT <- read.table("wt122data.tsv", header = TRUE, row.names = 1)
### Following replicates of the WT data were not used (Suzuki et al., 2018, BMC Genomics)
res <- c("04his3-12", "04his3-13", "04his3-18", "04his3-22", "04his3-33", "04his3-54", 
         "04his3-55", "04his3-60", "04his3-70", "04his385", "his3cnt-1", "his3old-1", "his3old-5")
AllHap$WT <- AllHap$WT[!rownames(AllHap$WT) %in% res,]
rm(res)
### Final data
AllHap$Data <- rbind(AllHap$WT, AllHap$Mut)


# Decoupling of non-linearity for the noise parameters
idx <- colnames(AllHap$Data)[grep("CV", colnames(AllHap$Data))]
names(idx) <- idx
idx <- sub("CV", "", idx)

#A matrix to save the output
SmoothSpan <- matrix(data = NA, nrow = length(idx), ncol = 2, dimnames = list(names(idx), c("Smooth span", "AIC")))

for(i in names(idx)){
  temp <- NULL
  aics <- NULL
  # x = Mean parameters
  # y = CV parameters
  temp <- data.frame(x=AllHap$Data[,idx[i]], y=AllHap$Data[,i])
  # AIC (test)
  aics <- sapply(10:99/100, function(i) AIC(gamlss(y~lo(~x, span=i), data=temp, family=NO)))
  SmoothSpan[i,"Smooth span"] <- (10:99/100)[which.min(aics)]
  SmoothSpan[i,"AIC"] <- min(aics)
}
rm(i)

save(SmoothSpan, file = "SmoothSpan_NonEssential.rdata")
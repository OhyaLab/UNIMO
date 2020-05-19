# Necessary code to calcualte the span for LOWESS regression

setwd("...")

set.seed(123)

# Prerequisite
if("gamlss" %in% rownames(installed.packages())){
  library(gamlss)
} else {install.packages("gamlss"); library(gamlss)}

Diploid_WT <- NULL
## Loading 114 diploid Wild-type strains:
### Average data: http://www.yeast.ib.k.u-tokyo.ac.jp/SCMD/download.php?path=wt114data.tsv
Diploid_WT$data <- read.table("wt114data.tsv", header = TRUE, row.names = 1)
colnames(Diploid_WT$data) <- chartr(".", "-", colnames(Diploid_WT$data))

# Decoupling of non-linearity for the noise parameters
idx <- colnames(Diploid_WT$data)[grep("CV", colnames(Diploid_WT$data))]
names(idx) <- idx
idx <- sub("CV", "", idx)

# A matrix to save the output
SmoothSpan <- matrix(data = NA, nrow = length(idx), ncol = 2, dimnames = list(names(idx), c("Smooth span", "AIC")))

for(i in names(idx)){
  temp <- NULL
  aics <- NULL
  # x = Mean parameters
  # y = CV parameters
  temp <- data.frame(x = Diploid_WT$data[,idx[i]], y = Diploid_WT$data[,i])
  # AIC (test)
  aics <- sapply(10:99/100, function(i) AIC(gamlss(y~lo(~x, span=i), data=na.omit(temp), family=NO)))
  SmoothSpan[i,"Smooth span"] <- (10:99/100)[which.min(aics)]
  SmoothSpan[i,"AIC"] <- min(aics)
}

save(SmoothSpan, file = "SmoothSpan_WT114.rdata")
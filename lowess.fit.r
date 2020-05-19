# args:
# x; input vector of x-axis (mean)
# y; input vector of y-axis (CV)
# ret; return type

lowess.fit <- function(x, y, ret="noise", f=0.2, fam="symmetric", deg=1, ite=4, sur="direct")
{
    ori <- data.frame(x=x, y=y)
    na.rm <- which(complete.cases(ori$x) & complete.cases(ori$y))
    d <- ori[na.rm,]
    fit <- loess(y~x, data=d, span=f, family=fam, degree=deg, control=loess.control(iterations=ite, surface=sur))  
    res <- NULL
    if(ret == "noise") {
        #browser()
        res0 <- ori$y[na.rm] - predict(fit)
        res <- rep(NA, nrow(ori))
        res[na.rm] <- res0
      } else if(ret == "fit") {
        res <- fit
      } else if(ret == "all") {
        res0 <- ori$y[na.rm] - predict(fit)
        res <- rep(NA, nrow(ori))
        res[na.rm] <- res0
        res <- list(fit=fit, res=res)
      }
    return(res)
}

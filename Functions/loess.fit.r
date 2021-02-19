# LOESS (locally estimated scatterplot smoothing) regression function
# to uncouple non-linear dependency between the coefficient of variation (CV) of their related mean (220 CalMorph parameters).
# Finally, instead of CV parameters, noise parameters as the residuals (i.e., observed value - predicted value) are calculated. 
## For a detailed description see Levy and Siegal (2008; PLoS Biol. 6(11):e264) and Yvert et al. (2013; BMC Syst. Biol. 7(1):54).

# Author: Shinsuke Ohnuki
##################################################
##################################################
##################################################
# Arguments are mainly explained in loess() function of 'stats' package [help(loess, package = "stats")]:
## x: Input vector of mean values
## y: Input vector of CV values
## ret: Return type
## f: Smooth span (to be defined in '01 NoiseParameters.r')
## fam: If "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function.
## deg: The degree of the polynomials to be used, normally 1 or 2.
## ite: The number of iterations used in robust fitting, i.e. only if 'fam' is "symmetric"
## sur: Should the fitted surface be computed exactly ("direct") or via interpolation from a kd tree ("interpolate")


loess.fit <- function(x, y, ret="noise", f=0.2, fam="symmetric", deg=1, ite=4, sur="direct")
{
    ori <- data.frame(x=x, y=y)
    na.rm <- which(complete.cases(ori$x) & complete.cases(ori$y))
    d <- ori[na.rm,]
    fit <- loess(y~x, data=d, span=f, family=fam, degree=deg, control=loess.control(iterations=ite, surface=sur))  
    res <- NULL
    if(ret == "noise") { # observed value - predicted value
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

# pBB() function of gamlss package (Stasinopoulos & Rigby, 2007; R. J. Stat. Softw., 23(7): 1-46)
# defines probability mass function of the beta-binomial distribution at any given sample
# between zero and one (0 < y < 1). In the case of very large quantiles (numerator equals denominator)
# very small p values (less than 1E-13) for the upper-tail test fail to be estimated.
# Following function (pBB2) solves this problem.

# Author Shinsuke Ohnuki
##################################################
##################################################
##################################################
# Arguments are explained in pBB() function of 'gamlss' package [help(pBB, package = "gamlss.dist")]:
## q: Vector of quantiles
## mu: Vector of positive probabilities
## Sigma: The dispersion parameter
## bd: Vector of binomial denominators
## lower.tail: Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
## log.p: Logical; if TRUE, probabilities p are given as log(p)

pBB2 <- function (q, mu = 0.5, sigma = 1, bd = 10, lower.tail = TRUE, 
                  log.p = FALSE) 
{
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1 ", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be >=0", "\n", ""))
  if (any(bd < q)) 
    stop(paste("y  must be <=  than the binomial denominator", 
               bd, "\n"))
  ly <- max(length(q), length(mu), length(sigma), length(bd))
  q <- rep(q, length = ly)
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)
  bd <- rep(bd, length = ly)
  if(lower.tail) {
    fn <- function(q, mu, sigma, bd) sum(dBB(0:q, mu = mu, sigma = sigma, 
                                             bd = bd))
  } else {
    fn <- function(q, mu, sigma, bd) sum(dBB(q:bd, mu = mu, sigma = sigma, 
                                             bd = bd))
  }
  Vcdf <- Vectorize(fn)
  cdf <- Vcdf(q = q, mu = mu, sigma = sigma, bd = bd)
  # Following three lines are omitted because it can produce P values less than zero,
  # when lower.tail equals FALSE and the numerator equals the denominator.
  # Instead of subtracting cdf from one in case of lower.tail=F, we added 36-42 lines to calculate
  # P values by accumulating the probability mass function as above.
  #
  # cdf <- if (lower.tail == TRUE)
  #     cdf
  # else 1 - cdf
  cdf <- if (log.p == FALSE) 
    cdf
  else log(cdf)
  if (length(sigma) > 1) 
    cdf2 <- ifelse(sigma > 1e-04, cdf, pBI(q, mu = mu, bd = bd, 
                                           lower.tail = lower.tail, log.p = log.p))
  else cdf2 <- if (sigma < 1e-04) 
    pBI(q, mu = mu, bd = bd, lower.tail = lower.tail, log.p = log.p)
  else cdf
  cdf2
}
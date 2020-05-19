# Solving the problem for very small p values for upper-tail test
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
# ifelse(q < bd, (q+1), q)
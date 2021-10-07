#' @title simulate normal distributed data
#' @description Simulate normal distributed data
#' @usage simulateData(pi0,m,n,rho, seed = NULL, power, alpha)
#' @param pi0 proportion of true null hypothesis
#' @param m number of variables, i.e., hypothesis
#' @param n number of observations
#' @param rho Level of equi-correlation between pairs of variables
#' @param seed specify seed 
#' @param power power
#' @param alpha alpha value for power analysis
#' @author Angela Andreella
#' @return Returns the simulated data matrix.
#' @export
#' @importFrom stats rnorm
#' @importFrom stats power.t.test

simulateData <- function(pi0,m,n, rho, seed = NULL, power = NULL, alpha = 0.05){
  if(is.null(seed)){set.seed(sample.int(1e5, 1))}
  if(is.null(n) & !is.null(power)){stop("Please insert sample size n")}
  m0 = round(m*pi0)
  m1 = round(m -m0)
  pwo <- power.t.test(power = power, n=n, sig.level = alpha, type = "one.sample", alternative = "two.sided", sd = 1)
  diff_mean<-pwo$delta
  #diff_mean <-rgamma(m1, shape = pwo$delta, scale = 1)
  #diff_mean <- rnorm(m1, mean = pwo$delta, sd = 1000)
  eps <- sqrt(1-rho)*matrix(rnorm(m*n), ncol=m)  + sqrt(rho)*matrix(rep(rnorm(n),m), ncol=m) 
  #mu <- c(diff_mean, rep(0, m0))
  mu <- c(rep(diff_mean,m1), rep(0, m0))
  X <- matrix(rep(mu, each=n), ncol=m) + eps
  
  return(t(X))
}



#' @title simulate normal distributed data
#' @description Simulate normal distributed data
#' @usage simulateData(pi0,m,rho,SNR, set.seed = NULL)
#' @param pi0 proportion of true null hypothesis
#' @param m number of variables, i.e., hypothesis
#' @param n number of observations
#' @param rho Level of equi-correlation between pairs of variables
#' @param set.seed specify seed 
#' @param d effect size
#' @param power power
#' @param alpha alpha value for power analysis
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Test} observed one sample t-test, \code{Test_H0} Test statistics under H0, \code{pv} observed p-values, \code{pv_H0} p-values under H0
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pwr.t.test

simulateData <- function(pi0,m,n, rho, set.seed = NULL, d = NULL, power = NULL, alpha = 0.1){
  if(is.null(set.seed)){set.seed(sample.int(1e5, 1))}
  if(is.null(d) & is.null(power)){stop("Please insert power or effect size")}
  if(is.null(n) & !is.null(power)){stop("Please insert sample size n")}
  m0 = round(m*pi0)
  m1 = round(m -m0)
  if(is.null(d)){d<-power.t.test(n = n, power = power, sig.level = alpha, type = "one.sample", alternative = "two.sided")$delta}
  sigma <- matrix(rep(rho,m*m),nrow = m,ncol=m) + diag(m)*(1-rho)
  eps <- rmvnorm(n = 1,mean = rep(0, nrow(sigma)),sigma = sigma)
  mu <- c(rep(0,m0),rep(d,m1))
  X <- mu + eps
  
  return(X)
}



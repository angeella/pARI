#' @title simulate normal distributed data
#' @description Simulate normal distributed data.
#' @usage simulateData(pi0,m,n, rho, seed = NULL, power = 0.8, alpha = 0.05)
#' @param pi0 Numeric value in `[0,1]`. Proportion of true null hypothesis.
#' @param m Numeric value. Number of variables. 
#' @param n Numeric value. Number of observations. 
#' @param rho Numeric value in `[0,1]`. Level of equi-correlation between pairs of variables.
#' @param seed Integer value. If you want to specify the seed. Default to to \code{NULL}
#' @param power Numeric value in `[0,1]`. Level of power. Default to 0.8.
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @author Angela Andreella
#' @return Returns a matrix with dimensions \eqn{m \times n}.
#' @export
#' @importFrom stats rnorm
#' @importFrom stats power.t.test

simulateData <- function(pi0,m,n, rho, seed = NULL, power = 0.8, alpha = 0.05){
  if(is.null(n) & !is.null(power)){stop("Please insert sample size n")}
  m0 = round(m*pi0)
  m1 = round(m -m0)
  if(is.null(seed)){set.seed(sample.int(1e5, 1))}else{set.seed(seed)}
  pwo <- power.t.test(power = power, n=n, sig.level = alpha, type = "one.sample", alternative = "two.sided", sd = 1)
  diff_mean<-pwo$delta
  #diff_mean <-rgamma(m1, shape = pwo$delta, scale = 1)
  #diff_mean <- rnorm(m1, mean = pwo$delta, sd = 1000)
  if(is.null(seed)){set.seed(sample.int(1e5, 1))}else{set.seed(seed)}
  eps <- sqrt(1-rho)*matrix(rnorm(m*n), ncol=m)  + sqrt(rho)*matrix(rep(rnorm(n),m), ncol=m) 
  #mu <- c(diff_mean, rep(0, m0))
  mu <- c(rep(diff_mean,m1), rep(0, m0))
  X <- matrix(rep(mu, each=n), ncol=m) + eps
  
  return(t(X))
}



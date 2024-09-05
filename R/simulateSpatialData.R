#' @title simulate normal distributed data 
#' @description Simulate normal distributed data with spatial correlation structure
#' @usage simulateSpatialData(pi0,m,n, theta, seed = NULL, power = 0.8, alpha = 0.05)
#' @param pi0 Numeric value in `[0,1]`. Proportion of true null hypothesis.
#' @param m Numeric value. Number of variables. 
#' @param n Numeric value. Number of observations. 
#' @param theta Numeric value in `[0,1]`. Level of correlation between pairs of variables. See details
#' @param seed Integer value. If you want to specify the seed. Default to to \code{NULL}
#' @param power Numeric value in `[0,1]`. Level of power. Default to 0.8.
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @author Angela Andreella
#' @return Returns a matrix with dimensions \eqn{m \times n}.
#' @description \code{theta} (\eqn{\theta}) describes how rapidly the correlation declines with respect to the distance between two voxels.
#' The three-dimensional coordinates of the voxels are defined as all combinations
#' of vector \eqn{c = (1, \dots, m1/3)}, then \eqn{\Sigma_\theta = \exp(-\theta K)} where \eqn{K} is the matrix containing the
#' euclidean distances between the three-dimensional coordinates' voxels.
#' So, \eqn{m^{1/3}} must be an integer value. 
#' @export
#' @importFrom stats rnorm
#' @importFrom stats power.t.test
#' @importFrom stats dist


simulateSpatialData <- function(pi0,m,n, theta, seed = NULL, power = 0.8, alpha = 0.05){
  if(is.null(n) & !is.null(power)){stop("Please insert sample size n")}
  m0 = round(m*pi0)
  m1 = round(m -m0)
  if(is.null(seed)){set.seed(sample.int(1e5, 1))}else{set.seed(seed)}
  pwo <- power.t.test(power = power, n=n, sig.level = alpha, type = "one.sample", alternative = "two.sided", sd = 1)
  diff_mean<-pwo$delta
  
  m3 <- round(m^(1/3),3)
  
  if(round(m3) != m3){stop("Please insert a valid number of variables")}
  
  simgrid_m <- expand.grid(1:m3, 1:m3, 1:m3)
  distance <- as.matrix(dist(simgrid_m))
  D <- chol(exp(-theta * distance))
  if(is.null(seed)){set.seed(sample.int(1e5, 1))}else{set.seed(seed)}
  eps <- matrix(rnorm(n * m), ncol = m) %*% D
  mu <- c(rep(diff_mean,m1), rep(0, m0))
  X <- matrix(rep(mu, each=n), ncol=m) + eps
  
  return(t(X))
}
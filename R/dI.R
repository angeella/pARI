#' @title Lower bound for the number of true discoveries
#' @description Calculates 1-\eqn{\alpha} lower confidence bound for the set-wise of false null hypotheses.
#' @usage dI(ix, cv, pvalues, iterative, approx, ncomb, ...)
#' @param ix Numeric vector: set-wise hypotheses considered. 
#' @param cv Numeric vector: critical vector computed by \code{\link{criticalVector}}.
#' @param pvalues If \code{iterative = TRUE} you must put here the matrix of \eqn{p}-values with 
#' dimensions \eqn{m \times B} where \eqn{m} is the number of variables and \eqn{B} the number of permutations. 
#' Instead, if \code{iterative = FALSE}, you can put directly the vector of \eqn{m} observed \eqn{p}-values.
#' @param iterative Boolean value. If \code{iterative = TRUE}, the iterative method is applied (computationally demanding). Default to \code{FALSE}. Please see the reference below.
#' @param approx Boolean value. Default to \code{TRUE}. If you are analyzing high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time. Please see the reference below.
#' @param ncomb Numeric value. If \code{approx = TRUE}, you must decide how many random sub collections (level of approximation) considered. Default to 100.
#' @param ... Further arguments for the iterative approach, i.e., \code{iterative = TRUE}.
#' @export
#' @author Angela Andreella
#' @return Numeric value: the lower confidence bound for the number of true discoveries concerning the cluster \code{ix} specified.
#' @importFrom utils combn
#' @references Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.
#' @examples
#'db <- simulateData(pi0 = 0.7, m = 100, n = 20, rho = 0)
#'out <- signTest(X = db)
#'pv <- cbind(out$pv, out$pv_H0)
#'cv <- criticalVector(pvalues = pv, family = "simes", lambda = 0.1, alpha = 0.1)
#'dI(ix = c(1:100), cv = cv, pvalues = pv)

dI <- function(ix, cv, pvalues, iterative = FALSE, approx = TRUE, ncomb = 100, ...){

  ## Default hidden arguments
  family=list(...)$family
  delta=list(...)$delta
  alpha=list(...)$alpha
  
  if(!iterative){
    family <- delta <- alpha <- NULL
  }
  if(iterative & !(exists("family", mode = "character") & exists("delta") & exists("alpha"))){
    stop("Please specify the family of confidence bounds, delta and alpha levels if you want to use the iterative approach")
  }
  if(!iterative){
    d <- permDiscoveries(ix = ix, cv = cv, praw = pvalues)
  }else{
    d <- permDiscoveries(ix = ix, cv = cv, praw = pvalues[,1])
  }
  
  if(iterative){
    d_seq <- c()
    d_seq[1] <- d
    it <- 1
    dist <- Inf
    while(dist !=0 & (d_seq[it] != length(ix))) {
      if(approx == TRUE){
        Kcomb <- replicate(ncomb, sample(ix,size = length(ix) - d_seq[it], replace=FALSE), 
                           simplify="matrix")
      }else{
        Kcomb <- combn(ix, length(ix) - d_seq[it]) 
        ncomb <- ncol(Kcomb)
      }
      #Create complementry set: combinations + all not in ix
      m <- dim(pvalues)[1]
      R <- which(!(c(1:m) %in% ix))
      lambda_kc <- sapply(c(1:ncomb), function(x) {
        if(is.matrix(Kcomb)){
          Kc <- Kcomb[,x]
        }else{Kc <- Kcomb[x]}
        Kc <- unique(c(Kc, R))
        P_Kc <- matrix(pvalues[Kc,], nrow= length(Kc), ncol = dim(pvalues)[2])
        lambdaCalibrate(X = P_Kc, family = family, alpha = alpha, delta = delta, m = dim(pvalues)[1])
      })
      lambda <- max(lambda_kc)
      cv <- criticalVector(pvalues= pvalues, family= family, 
                           alpha = alpha, delta = delta, lambda = lambda, m = dim(pvalues)[1])
      
      d_seq[it +1] <- permDiscoveries(ix = ix,cv = cv, praw = pvalues[,1])
      
      dist <- d_seq[it] - d_seq[it+1]
      it <- it + 1
      
    }
    d <- max(d_seq)
    
  }
  
  return(d)
}



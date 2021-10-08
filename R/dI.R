#' @title Lower bound for the number of true discoveries
#' @description Calculates (1-alpha) lower confidence bounds for the set-wise of false null hypotheses.
#' @usage dI(ix, cv, pvalues, iterative, approx, ncomb, family, alpha, delta)
#' @param ix numeric vector. It refers to the set-wise hypotheses considered. 
#' @param cv numeric vector. It refers to the critical vector computed by \code{\link{criticalVector}}.
#' @param pvalues matrix of pvalues with dimensions \eqn{m \times B} used instead of the data matrix \code{X}. Default to @NULL.
#' @param iterative Boolean value. If \code{iterative = TRUE}, the iterative method for improvement of confidence envelopes is applied. Default @FALSE.
#' @param approx Boolean value. Default @TRUE. If you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb Numeric value. If \code{approx = TRUE}, you must decide how many random subcollections (level of approximation) considered.
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector. 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param alpha numeric value in `[0,1]`. It expresses the alpha level to control the family-wise error rate.
#' @param delta numeric value. It expresses the delta value, please see the references. Default to 0. 
#' @export
#' @author Angela Andreella
#' @return numeric value: the lower confidence bound for the number of true discoveries concerning the cluster \code{ix} specified.
#' @importFrom utils combn
#' 
dI <- function(ix, cv, pvalues, iterative, approx, ncomb, family, alpha, delta){

  d <- permDiscoveries(ix = ix, cv = cv, praw = pvalues[,1])
  if(iterative){
    d_seq <- c()
    d_seq[1] <- d
    it <- 1
    dist <- Inf
    while(dist !=0) {
      if(approx == TRUE){
        Kcomb <- replicate(ncomb, sample(ix,size = length(ix) - d_seq[it], replace=FALSE), 
                           simplify="matrix")
      }else{
        Kcomb <- combn(ix, length(ix) - d_seq[it]) 
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



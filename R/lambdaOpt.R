#' @title Lambda calibration
#' @description \code{lambdaOpt} computes the optimal lambda calibration parameter used in the critical vector.
#' @usage lambdaOpt(pvalues, family, alpha, delta, step.down = FALSE, max.step = 10, m = NULL)
#' @param pvalues pvalues matrix with dimensions equal to the number of variables times the number of permutations.
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.#' @param alpha alpha level.
#' @param delta delta value. Do you want to consider sets with at least delta size? By default \code{delta = 0}. 
#' @param step.down by default \code{step.down = FALSE}. If you want to compute the lambda calibration parameter using the step down approach put \code{TRUE}.
#' @param max.step by default \code{max.step = 10}. Maximum number of steps for the step down approach.
#' @param m number of hypothesis, default is \code{NULL}.
#' @author Angela Andreella
#' @return lambda parameter estimate
#' @export
#' @importFrom stats pbeta

lambdaOpt <- function(pvalues, family, alpha, delta, step.down = FALSE, max.step = 10, m = NULL){
  #pvalues matrix with dimensions variables times permutations
  family_set <- c("simes", "aorc", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  if(is.null(m)){m <- dim(pvalues)[1]}

   ##TODO implementation step-down for beta
  lambdaE <- lambdaCalibrate(X = pvalues, alpha = alpha, delta = delta, family = family, m = m)
  

  if(step.down){
    convergence <- FALSE
    cv0 <- criticalVector(pvalues = pvalues, family = family, alpha = alpha, lambda = lambdaE, delta = delta)
    no_rej <- which(pvalues[,1] >= min(cv0))
    it <- 1
    lambda <- c()
    while(!convergence && it < max.step){
      lambda[it] <- lambdaCalibrate(X = pvalues[no_rej,], alpha = alpha, delta = delta, family = family, m = dim(pvalues)[1])
      cv1 <- criticalVector(pvalues = pvalues[no_rej,], family = family, alpha = alpha, lambda = lambda[it], delta = delta, m = dim(pvalues)[1])
      no_rej_new <- which(pvalues[,1] >= min(cv1))
      
      if(all(no_rej_new %in% no_rej)){
        convergence <- TRUE
        lambdaE <- lambda[it]
      }else{
        no_rej <- no_rej_new
      }
    }
  }

  
  return(lambdaE)
}



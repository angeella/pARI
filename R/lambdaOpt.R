#' @title Lambda calibration
#' @description \code{lambdaOpt} computes the optimal lambda calibration parameter used in the critical vector.
#' @usage lambdaOpt(pvalues, family, alpha = 0.05, delta = 0, step.down = FALSE,
#'  max.step = 10, m = NULL)
#' @param pvalues matrix of pvalues with dimensions \eqn{m \times B} used instead of the data matrix \code{X}. Default to @NULL.
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param alpha numeric value in `[0,1]`. It expresses the alpha level to control the family-wise error rate. Default 0.05.
#' @param delta numeric value. It expresses the delta value, please see the references. Default to 0. 
#' @param step.down Boolean value. Default @FALSE If you want to compute the lambda calibration parameter using the step-down approach put \code{TRUE}.
#' @param max.step Numeric value. Default to 10. Maximum number of steps for the step down approach, so useful when \code{step.down = TRUE}.
#' @param m numeric value. Number of hypothesis. Default @NULL.
#' @author Angela Andreella
#' @return numeric value. It expresses the lambda parameter estimate, plese see package references.
#' @export
#' @importFrom stats pbeta

lambdaOpt <- function(pvalues, family, alpha = 0.05, delta = 0, step.down = FALSE, max.step = 10, m = NULL){
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



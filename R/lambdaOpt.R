#' @title Lambda calibration
#' @description Computes the optimal lambda calibration parameter used in the critical vector \code{\link{criticalVector}}.
#' @usage lambdaOpt(pvalues, family, alpha = 0.05, delta = 0, step.down = FALSE,
#'  max.step = 10, m = NULL)
#' @param pvalues Matrix of \eqn{p}-values with dimensions \eqn{m \times B} where \eqn{m} is the number of variables 
#' and \eqn{B} the number of permutations used instead of the data matrix \code{X}. Default to \code{NULL}.
#' @param family String character. Name of the family confidence envelope to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"}, \code{"higher.criticism"}, and \code{"power"}.
#' Default to "simes".
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @param delta Numeric value. \eqn{\delta} value. Please see the reference below. Default to 0. 
#' @param step.down Boolean value. Default to \code{FALSE} If you want to compute the lambda calibration parameter using the step-down approach put \code{TRUE}. Please see the reference below.
#' @param max.step Numeric value. Default to 10. Maximum number of steps for the step down approach, so useful when \code{step.down = TRUE}.
#' @param m Numeric value. Number of hypotheses. Default to \code{NULL}.
#' @seealso \code{\link{criticalVector}}
#' @author Angela Andreella
#' @return Numeric value. \eqn{\lambda} parameter estimate.
#' @export
#' @references Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.
#' @importFrom stats pbeta
#' @examples 
#'db <- simulateData(pi0 = 0.8, m = 100, n = 20, rho = 0)
#'out <- signTest(X = db)
#'pv <- cbind(out$pv, out$pv_H0)
#'cv <- lambdaOpt(pvalues = pv, family = "simes", alpha = 0.05)


lambdaOpt <- function(pvalues, family = "simes", alpha = 0.05, delta = 0, step.down = FALSE, max.step = 10, m = NULL){
  #pvalues matrix with dimensions variables times permutations
  family_set <- c("simes", "aorc", "beta", "higher.criticism", "power")
  
  family <- match.arg(tolower(family), family_set)
  if(is.null(m)){m <- dim(pvalues)[1]}

  if(family == "beta"){
    pvalues[pvalues==0] <- .Machine$double.xmin
    
  }
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



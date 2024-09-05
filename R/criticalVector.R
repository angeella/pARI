#' @title Critical vector
#' @description Compute critical vector curve. 
#' @usage criticalVector(pvalues, family, alpha = 0.05, lambda, delta = 0, m = NULL)
#' @param pvalues Matrix of pvalues with dimensions \eqn{m \times B} used instead of the data matrix \code{X}. Default to @NULL.
#' @param family String character. Name of the family confidence envelope to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"}, \code{"higher.criticism"}, and \code{"power"}.
#' Default to "simes".
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @param lambda Numeric value. \eqn{\lambda} value computed by \code{\link{lambdaOpt}}.
#' Default to 1.
#' @param delta Numeric value. \eqn{\delta} value. Please see the reference below. Default to 0. 
#' @param m Numeric value. Number of hypothesis. Default to \code{NULL}.
#' @author Angela Andreella
#' @return Numeric vector. Critical vector curve with length \eqn{m}.
#' @export
#' @seealso \code{\link{lambdaOpt}}
#' @importFrom stats qbeta
#' @references Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.
#' @examples 
#'db <- simulateData(pi0 = 0.8, m = 100, n = 20, rho = 0)
#'out <- signTest(X = db)
#'pv <- cbind(out$pv, out$pv_H0)
#'cv <- criticalVector(pvalues = pv, family = "simes", lambda = 0.05)
#'plot(sort(pv[,1]), type = "l")
#'lines(cv)

criticalVector <- function(pvalues, family = "simes", alpha = 0.05, lambda, delta = 1, m = NULL){
  
  family_set <- c("simes", "aorc", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  #w <- dim(pvalues)[1]
  if(is.null(m)){m <- dim(pvalues)[1]}
  if(is.null(delta) ){delta = 0}
  if(family=="simes"){
    
    cv <- sapply(c(1:m), function(x) (((x-delta) * alpha * lambda)/(m-delta)))
  }
  if(family=="beta"){
    cv <- qbeta(lambda, c(1:m),m+1-c(1:m))
    cv <- unlist(sapply(c(1:length(cv)), function(x) if(is.na(cv[x])){qbeta(0, x,m+1-x)}else{cv[x]}))
  }
  if(family=="aorc"){
    
    cv <- sapply(c(1:m), function(x) (((x-delta) * lambda * alpha)/((m) - (x-delta) *(1 - lambda* alpha))))
    
  }

  if(family=="higher.criticism"){
    cv <- sapply(c(1:m), function(x) (2*x + lambda^2 - sqrt((2*x + lambda^2)^2 - 4*x^2 * (m + lambda^2)/m))/(2*(m + lambda^2))) 
  }
  if(family == "power"){
    cv <- sapply(c(1:m), function(x) (x/(m + sqrt(m)))^(-lambda))
  }

  return(cv)
}
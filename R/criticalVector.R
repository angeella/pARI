#' @title Critical value
#' @description compute critical values curve 
#' @usage criticalVector(pvalues, family, alpha, lambda, delta = NULL)
#' @param pvalues pvalues matrix with dimensiona variables times permutations
#' @param family Choose a family of confidence envelopes to compute the critical vector from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param alpha alpha
#' @param lambda lambda
#' @param delta delta
#' @author Angela Andreella
#' @return critical value curve
#' @export
#' @importFrom stats qbeta


criticalVector <- function(pvalues, family, alpha, lambda, delta = NULL){
  family_set <- c("simes", "aorc", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  #w <- dim(pvalues)[1]
  m <- dim(pvalues)[1]
  if(is.null(delta) ){delta = 0}
  if(family=="simes"){
    
    cv <- sapply(c(1:m), function(x) (((x-delta) * alpha * lambda)/(m-delta)))
  }
  if(family=="beta"){
    cv <- qbeta(lambda, c(1:m),m+1-c(1:m))
  }
  if(family=="aorc"){
    
    cv <- sapply(c(1:m), function(x) (((x-delta) * lambda * alpha)/((m) - (x-delta) *(1 - lambda* alpha))))
    
  }

  if(family=="higher.criticism"){
    #cv <- sapply(c(1:m), function(x) (2*x + lambda^2 - sqrt((2*x + lambda^2)^2 - 4*x^2 * (m + lambda^2)/m))/(2*(m + lambda^2))) 
    cv <- sapply(c(1:m), function(x) (2*x + lambda^2 - sqrt((2*x + lambda^2)^2 - 4*x^2 * (m + lambda^2)/m))/(2*(m + lambda^2))) 
    
  }
  
  return(cv)
}
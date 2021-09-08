#' @title Critical value
#' @description compute critical values curve 
#' @usage criticalVector(pvalues, family, alpha, lambda, delta = NULL)
#' @param pvalues pvalues matrix with dimensiona variables times permutations
#' @param family Choose a family of confidence envelopes to compute the critical vector from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param alpha alpha
#' @param lambda lambda
#' @param delta delta
#' @param m number of hypothesis
#' @author Angela Andreella
#' @return critical value curve
#' @export
#' @importFrom stats qbeta


criticalVector <- function(pvalues, family, alpha, lambda, delta = NULL, m = NULL){
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
    #cv <- sapply(c(1:m), function(x) (2*x + lambda^2 - sqrt((2*x + lambda^2)^2 - 4*x^2 * (m + lambda^2)/m))/(2*(m + lambda^2))) 
    #cv <- sapply(c(1:m), function(x) (2*x + lambda^2 - sqrt((2*x + lambda^2)^2 - 4*x^2 * (m + lambda^2)/m))/(2*(m + lambda^2))) 
    cv <- sapply(c(1:m), function(x) (2*x + lambda^2 + (lambda*sqrt(lambda^2*m + 4*m*x - 4*x^2)/sqrt(m)))/(2*(lambda^2 +m)))
  }
  
  return(cv)
}
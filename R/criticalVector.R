#' @title Critical vector
#' @description compute critical vector curve 
#' @usage criticalVector(pvalues, family, alpha, lambda, delta = NULL, m = NULL)
#' @param pvalues pvalues matrix with dimensions variables times permutations
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param alpha alpha level.
#' @param lambda lambda value computed by \code{\link{lambdaOpt}}.
#' @param delta delta value. Do you want to consider sets with at least delta size? By default \code{delta = 0}. 
#' @param m number of hypothesis, default is \code{NULL}.
#' @author Angela Andreella
#' @return critical vector curve
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
    cv <- sapply(c(1:m), function(x) (2*x + lambda^2 - sqrt((2*x + lambda^2)^2 - 4*x^2 * (m + lambda^2)/m))/(2*(m + lambda^2))) 
  }
  
  return(cv)
}
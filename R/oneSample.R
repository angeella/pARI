#' @title One-Sample t-test
#' @description Performs parametric one-sample t-tests.
#' @usage oneSample(X,alternative= "two.sided")
#' @param X data matrix where rows represent the \code{m} variables and columns the \code{n} observations.
#' @param alternative character string. It refers to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @author Angela Andreella
#' @return Returns a list with the following objects: 
#' - \code{Test}: vector with length equals \eqn{m}. Observed two-sample t-tests, one for each \eqn{m} variable, 
#' - \code{pv}: vector with length equals \eqn{m}. observed p-values, one for each \eqn{m} variable,
#' @export


oneSample<- function(X,alternative = "two.sided"){
  alternative_set <- c("two.sided", "greater", "lower")
  n <- ncol(X)
  m <- nrow(X)
  
  alternative <- match.arg(tolower(alternative), alternative_set)
  rowV <- rowVariance(X)

  Test <- ifelse(rowV==0,0, rowMeans(X)/(sqrt((rowV)/n)))
  pv <- switch(alternative, 
              "two.sided" = 2*(pt(abs(Test), df = n-1, lower.tail=FALSE)),
              "greater" = pt(Test, df = n-1, lower.tail=FALSE),
              "lower" = 1-pt(Test, df = n-1, lower.tail=FALSE))
  
  res <- list(Test = Test, pv = pv)
  
  return(res)
}
  

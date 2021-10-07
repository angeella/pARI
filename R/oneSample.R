#' @title One-Sample t-test
#' @description Performs parametric one-sample t-tests.
#' @usage oneSample(X,alternative= "two.sided")
#' @param X data matrix where rows represent the \code{m} variables and columns the \code{n} observations.
#' @param alternative character string referring to the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}). 
#' @author Angela Andreella
#' @return Returns a list with the following objects: Test test statistic, and pv corresponding raw pvalues.
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
  

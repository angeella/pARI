#' @title One-Sample t-test
#' @description performs the One-Sample t-test for each feature.
#' @usage oneSample(X,alternative= "two.sided")
#' @param X data matrix where rows represents the variables and columns the observations.
#' @param alternative a character string referring to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @author Angela Andreella
#' @return Returns a list with the following objects: Test test statistic, and pv corresponding raw pvalues.
#' @export


oneSample<- function(X,alternative = "two.sided"){
  alternative_set <- c("two.sided", "greater", "lower")
  n <- ncol(X)
  m <- nrow(X)
  
  alternative <- match.arg(tolower(alternative), alternative_set)
  rowV <- rowVariance(X)
  #rowV <- ifelse(rowV==0,.Machine$double.xmin, rowV)
  #Test = rowMeans(X)/(sqrt((rowV)/n))
  Test <- ifelse(rowV==0,0, rowMeans(X)/(sqrt((rowV)/n)))
  pv <- switch(alternative, 
              "two.sided" = 2*(pt(abs(Test), df = n-1, lower.tail=FALSE)),
              "greater" = pt(Test, df = n-1, lower.tail=FALSE),
              "lower" = 1-pt(Test, df = n-1, lower.tail=FALSE))
  
  res <- list(Test = Test, pv = pv)
  
  return(res)
}
  

#' @title One Sample t test
#' @description Performs One Sample t test
#' @usage oneSample(X,alternative)
#' @param X data where rows represents the variables and columns the observations
#' @param alternative = character referring to the alternative hypothesis, "two.sided", "greater" or "less". Default is "two.sided"
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Test} test statistic, and \code{pv} corresponding raw pvalues.
#' @export



oneSample<- function(X,alternative = "two.sided"){
  alternative_set <- c("two.sided", "greater", "lower")
  
  n <- ncol(X)
  m <- nrow(X)
  
  alternative <- match.arg(tolower(alternative), alternative_set)
  
  Test = rowMeans(X)/(sqrt((rowVariance(X))/n))
  pv <- switch(alternative, 
              "two.sided" = 2*(pnorm(abs(Test), lower.tail=FALSE)),
              "greater" = pnorm(Test, lower.tail=FALSE),
              "less" = 1-pnorm(Test, lower.tail=FALSE))
  
  res <- list(Test = Test, pv = pv)
  
  return(res)
}
  

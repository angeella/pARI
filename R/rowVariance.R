#' @title Rows Variance
#' @description performs the variance for each row in a matrix. Internal function.
#' @usage rowVariance(X,na.rm = TRUE) 
#' @param X data matrix where rows represent the \eqn{m} variables and columns the \eqn{n} observations.
#' @param na.rm Boolean value. Default @TRUE. 
#' @author Angela Andreella
#' @return vector with length equals \eqn{m} rows variance.

rowVariance <- function (X,na.rm = TRUE) 
{
  sqr = function(X) X * X
  n = rowSums(!is.na(X))
  n[n <= 1] = NA
  return(rowSums(sqr(X - rowMeans(X,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}
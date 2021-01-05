#' @title Rows Variance
#' @description performs the variance for each row in a matrix.
#' @usage rowVariance(X,na.rm = TRUE) 
#' @param X data matrix where rows represents the variables and columns the observations.
#' @param na.rm by default \code{na.rm=TRUE}. 
#' @author Angela Andreella
#' @return rows variance
#' @export

rowVariance <- function (X,na.rm = TRUE) 
{
  sqr = function(X) X * X
  n = rowSums(!is.na(X))
  n[n <= 1] = NA
  return(rowSums(sqr(X - rowMeans(X,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}
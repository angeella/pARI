#' @title Rows Variance
#' @description Performs the variance for each row in a matrix
#' @param X data where rows represents the variables and columns the observations
#' @param na.rm remove na? 
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
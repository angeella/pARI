#' @title tdp
#' @description compute tdp
#' @param hommel hommel
#' @param ix ix
#' @param alpha alpha
#' @author hommel package
#' @return tdp
#' @export

tdp <- function (hommel, ix, alpha) 
{
  m <- length(hommel@p)
  if (missing(ix)) {
    d <- discoveries(hommel, alpha = alpha)
    k <- m
  }
  else {
    p <- hommel@p[ix]
    k <- length(p)
    d <- discoveries(hommel, ix, alpha = alpha)
  }
  d/k
}
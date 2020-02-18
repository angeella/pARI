#' @title Higher Criticism 
#' @description Performs HC by Donoho and Jim
#' @param p = ordered raw pvalues




HC <- function(p){
  
  m <- length(p)
  i <- c(1:m)
  hc <- sqrt(m) * (i/m - p)/(p*(1-p))
  
  return(hc)
}
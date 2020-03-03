#' @title True Discovery Permutation
#' @description performs lower confidence bound of true discoveries using permutation local test.
#' @param ix The selection of hypotheses considered
#' @param cv critical values 
#'  @param p raw pvalues

TDPerm <- function(ix,cv,p){
  discoveries <- dI(ix,cv,p)
  return(list(discoveries= discoveries))
}
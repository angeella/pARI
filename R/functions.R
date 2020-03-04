#' @title n rejections using permutation theory
#' @description Compute the number of rejections having pvalues computed using permutation distribution
#' @usage rejPerm(pvalues,ct)
#' @param pvalues matrix with dimension number of permutation times number of tests
#' @param ct set of cutoff between 0 and 1
#' @author Angela Andreella
#' @return Returns number of rejections
#' @export

rejPerm <- function(pvalues,ct){
  rejs <- apply(pvalues,1,function(i) sum(i < ct))
  return(rejs)
}


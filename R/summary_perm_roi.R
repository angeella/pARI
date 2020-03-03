#' @title Summary permutation ROI
#' @description Performs Single step closed testing using permutation local test
#' @param cv = vector of critical values
#' @param ix = vector of index regarding the set of hypothesis to analyze
#' @param pvalues = raw pvalues
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Size} Total number of hypotheses rejected, \code{FalseNull} number of false null hypotheses, \code{TrueNull} number of true null hypotheses, and \code{ActiveProp} active proportion inside ix.
#' @export



summary_perm_roi <- function(cv,ix,pvalues){
  idix <- which(ix)
  p <- pvalues[1,idix]
  Total = length(p)
  False_Null= dI(ix = idix,cv = cv,praw = pvalues[1,])
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  list(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
}

# @description computes alpha-level estimate of significant variables
# @pvalues = pvalues matrix

summary_perm_roi <- function(cv,ix,pvalues, iterative, approx, ncomb, family, delta, alpha){
  idix <- which(ix)
  p <- pvalues[idix,1]
  Total = length(idix)
  False_Null= dI(ix = idix,cv = cv,pvalues = pvalues, iterative, approx, ncomb, family, delta, alpha)
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  list(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
}

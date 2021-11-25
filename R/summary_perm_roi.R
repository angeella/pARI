# @description Internal function.
# @pvalues = pvalues matrix

summary_perm_roi <- function(cv,ix,pvalues, iterative, approx, ncomb, ...){
  idix <- which(ix)
  p <- pvalues[idix,1]
  Total = length(idix)
  False_Null= dI(ix = idix,cv = cv,pvalues = pvalues, 
                 iterative = iterative, approx = approx, ncomb = ncomb, ...)
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  list(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
}

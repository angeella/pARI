# @description computes alpha-level estimate of significant variables
# @pvalues = raw pvalues

summary_perm_roi <- function(cv,ix,pvalues){
  idix <- which(ix)
  p <- pvalues[idix]
  Total = length(idix)
  False_Null= dI(ix = idix,cv = cv,praw = pvalues)
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  list(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
}

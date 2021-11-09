# @description Internal function.
# @p = raw p-values

TDPerm <- function(ix,cv,p){
  
  discoveries <- dI(ix = ix,cv = cv,pvalues = p)
  
  return(list(discoveries= discoveries))
}
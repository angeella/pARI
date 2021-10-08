# @description Internal function.
# @p = raw p-values

TDPerm <- function(ix,cv,p){
  
  discoveries <- dI(ix,cv,p)
  
  return(list(discoveries= discoveries))
}
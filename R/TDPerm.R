# @description computes the lower bound for the number of true discoveries inside the set of variables selected.
# @p = raw p-values

TDPerm <- function(ix,cv,p){
  
  discoveries <- dI(ix,cv,p)
  
  return(list(discoveries= discoveries))
}
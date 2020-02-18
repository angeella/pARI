#True Discovery Proportion Permutation
##
# @param ix The selection of hypotheses considered
# @param cv critical values 
# @param p raw pvalues

TDPerm <- function(ix,cv,p){
  discoveries <- dI(ix,cv,p)
  return(list(discoveries= discoveries))
}
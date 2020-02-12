############################FUNCTIONS#######################

#NUMBER OF REJECTIONS UNDER PERMUTATION

#pvalues= matrix with dimension w (number of permutation) times m (number of tests)
#ct = cutoff
rejPerm <- function(pvalues,ct){
  rejs <- apply(pvalues,1,function(i) sum(i < ct))
  return(rejs)
}


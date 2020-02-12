#######################Critical values #########################

cv <- function(pvalues, family, alpha, shift = NULL, lambda, ct = NULL, delta = NULL){
  
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  if(family=="Simes"){
    if(is.null(shift) ){shift = 0}
    cv <- sapply(c(1:m), function(x) ((x * alpha * lambda)/m)- shift)
  }
  if(family=="Beta"){
    cv <- qbeta(lambda, c(1:m),m+1-c(1:m))
  }
  if(family=="Finner"){
    
    cv <- sapply(c(1:m), function(x) ((x * lambda * alpha)/(m - x *(1 - lambda* alpha))))
    
  }
  #if(family=="SimesCluster"){
  #  cv <- sapply(c(1:m), function(x) (((x-8) * alpha * lambda)/m))
  #  
  #  }
  if(family=="SimesCluster"){
    #cl <- max(sapply(c(1:w), function(x) sum(pvalues[x,] <= min(ct)))) +8
    #cv <- sapply(c(1:m), function(x) (((x-400) * alpha * lambda)/m))
    cv <- sapply(c(1:m), function(x) (((x-delta) * alpha * lambda)/(m-delta)))
  }
  
  if(family=="FinnerCluster"){
    #cl <- max(sapply(c(1:w), function(x) sum(pvalues[x,] <= min(ct)))) +8
    #cv <- sapply(c(1:m), function(x) (((x-400) * alpha * lambda)/m))
    cv <- sapply(c(1:m), function(x) (((x-delta) * lambda * alpha)/(m - (x-delta) *(1 - lambda* alpha))))
  } 
  
  
  return(cv)
}
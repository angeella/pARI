#######################Critical values #########################

cv <- function(pvalues, family, alpha, lambda, ct = NULL, delta = NULL){
  
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  if(is.null(delta) ){delta = 0}
  if(family=="Simes"){
    
    cv <- sapply(c(1:m), function(x) (((x-delta) * alpha * lambda)/(m-delta)))
  }
  if(family=="Beta"){
    cv <- qbeta(lambda, c(1:m),m+1-c(1:m))
  }
  if(family=="Finner"){
    
    cv <- sapply(c(1:m), function(x) (((x-delta) * lambda * alpha)/(m - (x-delta) *(1 - lambda* alpha))))
    
  }

  
  
  return(cv)
}
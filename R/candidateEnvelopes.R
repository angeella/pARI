#########################FAMILIES OF CANDIDATE ENVELOPES##################################

############SIMES############


simesB <- function(cutoff, lambda, m, delta, alpha){
  quants <- lapply(c(1:m), function(i) {lambda*i*alpha/m -delta} )
  out <- sum(quants <= cutoff)
  return(out)
}

############BETA############

betaB <- function(cutoff, lambda, m){
  quants <- lapply(c(1:m), function(i) {qbeta(lambda,i,m+1-i)} )
  out <- sum(quants <= cutoff)
  return(out)
  
}

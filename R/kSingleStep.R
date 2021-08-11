kSingleStep <- function(pvalues, family, alpha, delta, ix, iterative = FALSE, approx = TRUE, step.down = FALSE, max.step = 100){
 # source("R/lambdaR.R")
  #source("src/rowSortC.cpp")
  
  lambda <- lambdaOpt(pvalues = pvalues, family = family, alpha = alpha, delta = delta)
  cv <- criticalVector(pvalues=pvalues, family= family, alpha = alpha, delta = delta, lambda = lambda)
  discoveries <- c()
  discoveries[1] <- dI(ix =  ix, cv = cv, pvalues = pvalues, iterative = FALSE)
  M <- nrow(pvalues)
  it <- 1
  dist <- Inf
  while(dist <= 0){
    it <- it + 1
    if(it ==2){
      k <- discoveries[it-1]
    }else{
      k <- discoveries[it-1] - discoveries[it-2]
      
    }
    m <- length(ix)
    ix_c <- which(!(c(1:M) %in% ix))
    
 #   worst_p <- sort(pvalues[ix,1])[c((k+1):m)]
    worst_p_idx <- order(pvalues[ix,1])[c((k+1):m)]
    worst_p_idx <- ix[worst_p_idx]
  #  best_p <- sapply(c(2:ncol(pvalues)), function(x) sort(pvalues[ix,x])[c(1:(m-k))])
    best_p_idx <- sapply(c(2:ncol(pvalues)), function(x) c(ix[order(pvalues[ix,x])[c(1:(m-k))]], ix_c))
    
    P <- cbind(P[,1][c(worst_p_idx, ix_c)], 
               matrix(P[,2:ncol(pvalues)][c(best_p_idx)], ncol = ncol(pvalues)-1, nrow = M-k))
    lambda <- lambdaOpt(pvalues = P, family = family, alpha = alpha, delta = delta)
    cv <- criticalVector(pvalues=pvalues, family= family, alpha = alpha, delta = delta, lambda = lambda)
    discoveries[it] <- dI(ix=ix,cv = cv,pvalues = pvalues, iterative = FALSE) 
    dist <- discoveries[it] - discoveries[it-1]
  }
  
  print(discoveries)
  
  
  discoveries = discoveries[it]
  
  return(discoveries)
}
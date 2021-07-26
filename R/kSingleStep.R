kSingleStep <- function(P, family, alpha, delta){
  source("R/lambdaR.R")
  source("src/rowSortC.cpp")
  
  lambda <- lambdaOpt(P, family = family, alpha = alpha, delta = delta)
  cv <- criticalVector(pvalues=P, family= family, alpha = alpha, delta = delta, lambda = lambda)
  discoveries <- c()
  discoveries[1] <- dI(c(1:dim(P)[1]),cv,P, iterative = FALSE)
  it <- 1
  dist <- Inf
  while(dist != 0){
    it <- it + 1
    if(it ==2){
      k <- discoveries[it-1]
    }else{
      k <- discoveries[it-1] - discoveries[it-2]
      
    }
    m <- nrow(P)
    worst_p <- sort(P[,1])[c((k+1):m)]
    worst_p_idx <- order(P[,1])[c((k+1):m)]
    best_p <- sapply(c(2:ncol(P)), function(x) sort(P[,x])[c(1:(m-k))])
    best_p_idx <- sapply(c(2:ncol(P)), function(x) order(P[,x])[c(1:(m-k))])
    
    P <- cbind(P[,1][worst_p_idx], matrix(P[,2:ncol(P)][best_p_idx], ncol = ncol(P)-1, nrow = m-k))
    lambda <- lambdaOptR(P, family = family, alpha = alpha, delta = delta)
    cv <- criticalVector(pvalues=P, family= family, alpha = alpha, delta = delta, lambda = lambda)
    discoveries[it] <- max(sapply(c(1:length(cv)), function(x) 1 - x + sum(P[,1] <= cv[x]))) + discoveries[it-1]
    dist <- discoveries[it] - discoveries[it-1]
    print(discoveries[it])
  }
  
  
  
  discoveries = discoveries[it]
  
  return(discoveries)
}
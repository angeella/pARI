
permDiscoveriesIt <- function(ix, cv, pvalues, approx = TRUE, ncomb, family, alpha, delta){
  
 # family_set <- c("simes", "aorc", "beta", "higher.criticism")
  #family <- match.arg(tolower(family), family_set)
  set.seed(NULL)
  
  #Single-step
  Bt <- c()
  m <- dim(pvalues)[1]
 lambda <- lambdaCalibrate(X = pvalues, family = family, alpha = alpha, delta = delta)
  #lambda <- lambdaOptR(pvalues = pvalues, family = family, alpha = alpha, delta = delta)
  cv <- criticalVector(pvalues = pvalues, family= family, alpha = alpha, delta = delta, lambda = lambda)
 # print("1")
 # print(dim(pvalues))
  Bt[1] <- length(ix) - permDiscoveries(ix = ix,cv = cv,praw = pvalues[,1])
  
  if(Bt[1] == 0){
    #If I don t have False discoveries
    discoveries <- length(ix)
  }else{
      it <- 1
      dist <- Inf
      while(dist !=0) {
        if(approx == TRUE){
          Kcomb <- replicate(ncomb, sample(ix,size =Bt[it], replace=FALSE), simplify="matrix")
        }else{
          Kcomb <- combn(ix, Bt[it]) 
        }
        #Create complementry set: combinations + all not in ix
        
        R <- which(!(c(1:m) %in% ix))
        print("2")
        lambda_kc <- sapply(c(1:ncomb), function(x) {
          if(is.matrix(Kcomb)){
            Kc <- Kcomb[,x]
          }else{Kc <- Kcomb[x]}
          Kc <- unique(c(Kc, R))
          P_Kc <- matrix(pvalues[Kc,], nrow= length(Kc), ncol = dim(pvalues)[2])
          lambdaOptR(pvalues = P_Kc, family = family, alpha = alpha, delta = delta)
        })
        print("3")
        lambda <- max(lambda_kc)
        cv <- criticalVector(pvalues= pvalues, family= family, alpha = alpha, delta = delta, lambda = lambda)
        Bt[it +1] <- length(ix) - permDiscoveries(ix = ix,cv = cv,praw = pvalues[,1])
        
        dist <- Bt[it] - Bt[it+1]
        it <- it + 1
        
      }
      print("4")
    print(Bt)
    B_est <- min(Bt)
    
    discoveries <- length(ix) - B_est
  }
  
  return(discoveries)
}
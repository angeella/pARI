# library(pARI)
# X <- simulateData(pi0 = 0.9, m = 20, n = 10, rho = 0, power = 0.9)
# ix <- c(1:19)
# alpha = 0.05
# family = "Simes"
# delta = 0
# B = 1000
# test.type = "one_sample"
# ncomb = 2
# approx = TRUE
# #X: dati variables times obs
# 
# 
# out <- pARIt(X = X, ix = ix, ncomb = 200, approx = TRUE)
# out
# 
# out <- pARIt(X = X, ix = ix, ncomb = 200, approx = FALSE)
# out

pARIt <- function(X, ix, alpha = 0.05, family = "Simes", 
                  delta = 0, B = 1000, test.type = "one_sample", approx = TRUE, ncomb, ...){
  

  Bt <- c()
  Bt[1] <- length(ix) - pARI(X = X,alpha = alpha, test.type = test.type,
                             family = family, ix = ix)$discoveries
  converge <- FALSE
  it <- 2
  while (converge == FALSE) {
    
    if(approx == TRUE){
    
      Kc1 <- replicate(ncomb, sample(ix,size =  length(ix) - Bt[it-1], replace=FALSE), simplify="matrix")
      if(sum(!(c(1:nrow(X)) %in% ix))!=0){
        Kc2 <- replicate(ncomb, sample(which(!(c(1:nrow(X)) %in% ix)),size =  length(ix) - Bt[it-1], replace=FALSE), simplify="matrix")
        Kcomb <- cbind(Kc1,Kc2)
        }else{
        Kcomb <- Kc1
      }
      
    }else{
    Kcomb <- combn(ix, length(ix) - Bt[it-1]) 
  
    }
    
    B_kc <- sapply(c(1:ncol(Kcomb)), function(x) 
      length(ix[(!(ix %in% Kcomb[,x]))]) - pARI(X = X,alpha = alpha, test.type = test.type,
                                                family = family, ix = ix[(!(ix %in% Kcomb[,x]))])$discoveries
    )
    
    Bt[it] <- max(B_kc)
    
    if(it!=2 & Bt[it] == Bt[it-1]){
      converge <- TRUE
    }
    it <- it+ 1
    #print(it)
  }
  
  B_est <- min(Bt)
  discoveries <- length(ix) - B_est
  return(discoveries)
}

permDiscoveriesIt <- function(ix, cv, praw, approx = TRUE, ncomb){
  
  
  Bt <- c()
  Bt[1] <- length(ix) - permDiscoveries(ix = ix,cv = cv,praw = praw)
  if(Bt[1] == 0){
    discoveries <- length(ix)
  }else{
    if(length(ix) - Bt[1] !=0){
      
      converge <- FALSE
      it <- 2
      while (converge == FALSE) {
        
        if(approx == TRUE){
          Kcomb <- replicate(ncomb, sample(ix,size =  length(ix) - Bt[it-1], replace=FALSE), simplify="matrix")
        }else{
          Kcomb <- combn(ix, length(ix) - Bt[it-1]) 
          
        }
        B_kc <- sapply(c(1:ncol(Kcomb)), function(x) 
          length(ix[(!(ix %in% Kcomb[,x]))]) - permDiscoveries(ix = ix[(!(ix %in% Kcomb[,x]))],cv = cv,praw = praw)
        )
        
        
        
        
        Bt[it] <- max(B_kc)
        print(Bt)
        if(it!=2 & Bt[it] == Bt[it-1]){
          converge <- TRUE
        }
        it <- it+ 1
        #print(it)
      }
    }
    

    
    B_est <- min(Bt)
    discoveries <- length(ix) - B_est
  }
  
  return(discoveries)
}
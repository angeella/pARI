
permDiscoveriesIt <- function(ix, cv, praw, approx = TRUE, ncomb){
  
  set.seed(NULL)
  Bt <- c()
  Bt[1] <- length(ix) -permDiscoveries(ix = ix,cv = cv,praw = praw)
  
  if(Bt[1] == 0){
    discoveries <- length(ix)
  }else{
    if(length(ix) - Bt[1] !=0){
      
      converge <- TRUE
      it <- 2
      while(converge) {
        
        if(approx == TRUE){
          Kcomb <- replicate(ncomb, sample(ix,size =  length(ix) - Bt[it-1], replace=FALSE), simplify="matrix")
        }else{
          Kcomb <- combn(ix, length(ix) - Bt[it-1]) 
        }
        # if(it == 2) {print(Kcomb)}
    
        nK <- ifelse(is.matrix(Kcomb), ncol(Kcomb), length(Kcomb))
        Bt_kc <- sapply(c(1:nK), function(x) {
          if(is.matrix(Kcomb)){
            Kc <- which(!(ix %in% Kcomb[,x]))
          }else{Kc <- which(!(ix %in% Kcomb[x]))}
          length(Kc) - permDiscoveries(ix = Kc,cv = cv,praw = praw)})
        
        
        Bt[it] <- max(Bt_kc)
        print(Bt[it])
        if(Bt[it] == Bt[it -1]){
          converge <- FALSE
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
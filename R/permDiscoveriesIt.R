
permDiscoveriesIt <- function(ix, cv, praw, approx = TRUE, ncomb){
  
  set.seed(NULL)
  
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
          Kcomb <- replicate(ncomb, sample(ix,size =  Bt[it-1], replace=FALSE), simplify="matrix")
        }else{
          Kcomb <- combn(ix, Bt[it-1]) 
          
        }
        # if(it == 2) {print(Kcomb)}
        
        Rc <- which(!(c(1:length(praw)) %in% ix))
        
        B_kc <- sapply(c(1:ncol(Kcomb)), function(x) {
         # Kc <- ix[!((ix %in% Kcomb[,x]))]
          Kc <- Kcomb[,x]
          SetC <- c(Rc, Kc)
          length(Kc) - permDiscoveries(ix = Kc,cv = cv,praw = praw)})
        
        
        Bt[it] <- max(B_kc)
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
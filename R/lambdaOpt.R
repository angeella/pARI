###########################Lambda Optimal############################

#For every permutation, I compute the best possible lambda. Then I take the alpha w +1quantile due to discreteness.

#The idea with the min and max cut-offs is that you only take into account the p-values between 
#these cut-offs. So the smaller the interval (mincutoff, maxcutoff) is, the more power you will have. 
#The interval (mincutoff, maxcutoff) corresponds with the set \Gamma in the paper.

#pvalues = ordered pvalues 
#ct= vector of cutoffs
#family = Beta or Simes or Finner

lambdaOpt <- function(pvalues, family, ct, alpha, delta = NULL){
  l <- c()
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  if(is.null(delta) ){delta = 0}
  for(j in 1:w){
    minc <- sum(pvalues[j,] <=min(ct)) + 1
    maxc <- sum(pvalues[j,] <=max(ct))
    if(family =="Simes"){
      
      minc = minc + delta
      lambda <- ((m-delta)*(pvalues[j,minc:maxc]))/((c(minc:maxc)-delta)*alpha)
      #lambda <- (m*(pvalues[j,minc:maxc] + shift))/(c(minc:maxc))
    }
    if(family == "Beta"){
      lambda <- pbeta(pvalues[j,c(minc:maxc)],c(minc:maxc),m+1-c(minc:maxc))
      
    }
    
    if(family == "Finner"){
      minc = minc + delta
      lambda <- (pvalues[j,c(minc:maxc)]*(m - c(minc:maxc) + delta) )/ (alpha * (c(minc:maxc) - delta - (pvalues[j,c(minc:maxc)]* c(minc:maxc) + (delta*pvalues[j,c(minc:maxc)]))))
      
    }

    l[j] <- min(lambda)
  }
  
  lambdaE <- sort(l)[floor(alpha*w)+1]
  
  #nRej_Perm <- rejPerm(pvalues,min(ct))
  #quantRej_min <- sort(nRej_Perm)[ceiling((1-alpha)*w)]
  #if(quantRej_min>0){
  #lambdaE <- min(lambdaE, (min(ct) + shift)*m/(quantRej_min*alpha))}
  
  
  
  return(lambdaE)
}


lambdaOptAprox <- function(pvalues, family, ct = tS, alpha, shift = delta,cb, Kc){
  l <- c()
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  quant <- (pvalues+shift)*m/(alpha)
  lambdaA <- c()
  for(j in 1:w){
    minc <- sum(pvalues[j,Kc[,cb]] <=min(ct)) + 1
    maxc <- sum(pvalues[j,Kc[,cb]] <=max(ct))
    
    if(is.null(shift)){shift = 0}
    lambdaA[j] <- min(sapply(c(minc:maxc), function(x) quant[j, Kc[   sort.list(pvalues[j,Kc[,cb]])[x] ,cb   ] ]/x  ))
    #lamb <- min(lamb,  betaquantsS[, combs2[   sort.list(pvmatr.uns[j,combs2[,c]]),c]][j,a]/a  )
    
    
  }
  
  
  # nRej_Perm <- rejPerm(pvalues,min(ct))
  # quantRej_min <- sort(nRej_Perm)[ceiling((1-alpha)*w)]
  # if(quantRej_min>0){
  #   lambdaE <- min(lambdaE, (min(ct) + shift)*m/(quantRej_min*alpha))}
  
  return(lambdaA)
}

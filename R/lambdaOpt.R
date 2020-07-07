#' @title Lambda calibration
#' @description compute lambda parameter
#' @usage lambdaOpt(pvalues,family,ct,alpha, delta = NULL)
#' @param pvalues pvalues matrix
#' @param family family
#' @param alpha alpha
#' @param ct set threshold
#' @param delta delta
#' @author Angela Andreella
#' @return lambda
#' @export
#' @importFrom stats pbeta


lambdaOpt <- function(pvalues, family, ct = c(0,1), alpha, delta = NULL){
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  if(is.unsorted(pvalues[1,])){pvalues = rowSortC(pvalues)}
  l <- c()
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  if(is.null(delta) ){delta = 0}
  for(j in 1:w){
    minc <- sum(pvalues[j,] <=min(ct)) + 1
    maxc <- sum(pvalues[j,] <=max(ct))
    if(family =="simes"){
      
      minc = minc + delta
      lambda <- ((m-delta)*(pvalues[j,minc:maxc]))/((c(minc:maxc)-delta)*alpha)
      #lambda <- (m*(pvalues[j,minc:maxc] + shift))/(c(minc:maxc))
    }
    if(family == "beta"){
      lambda <- pbeta(q =pvalues[j,c(minc:maxc)],shape1 = c(minc:maxc),shape2 =m+1-c(minc:maxc))
      
    }
    
    if(family == "finner"){
      minc = minc + delta
      #lambda <- (pvalues[j,c(minc:maxc)]*(m - c(minc:maxc) + delta) )/ (alpha * (c(minc:maxc) - delta - (pvalues[j,c(minc:maxc)]* (c(minc:maxc) - 1))))
      lambda <- (pvalues[j,c(minc:maxc)]*(m - 1) )/ (alpha * (c(minc:maxc) - delta) * (1 - pvalues[j,c(minc:maxc)]))
      
    }
    
    
    if(family =="higher.criticism"){
      
      lambda <- (sqrt(m)*((c(minc:maxc)/m) - pvalues[j,c(minc:maxc)]))/(sqrt(pvalues[j,c(minc:maxc)]*(1-pvalues[j,c(minc:maxc)])))
      
      #HigherCriticism <- function(lambda,pvalue){
        
      # m <- length(pvalue)
      #  i <- c(1:m)
      #  (2*i + lambda^2 - sqrt((2*i + lambda^2)^2 - 4*i^2 * (m + lambda^2)/m))/(2*(m+ lambda^2)) - pvalue
        
      #}
      #lambda <- sapply(c(minc:maxc), function(x) rootSolve::uniroot.all(HigherCriticism,pvalue = pvalues[j,x],lower=0, upper=1000))
      #lambda <- unlist(lambda)
}
    if(family=="beta"){
      l[j] <- min(lambda[lambda>0])
    }else{
      l[j] <- min(lambda)
    }
    
  }
  
  lambdaE <- sort(l)[floor(alpha*w)+1]

  #nRej_Perm <- rejPerm(pvalues,min(ct))
  #quantRej_min <- sort(nRej_Perm)[ceiling((1-alpha)*w)]
  #if(quantRej_min>0){
  #lambdaE <- min(lambdaE, (min(ct) + shift)*m/(quantRej_min*alpha))}
  
  
  
  return(lambdaE)
}


#lambdaOptAprox <- function(pvalues, family, ct = tS, alpha, shift = delta,cb, Kc){
#  l <- c()
#  w <- dim(pvalues)[1]
#  m <- dim(pvalues)[2]
#  quant <- (pvalues+shift)*m/(alpha)
#  lambdaA <- c()
#  for(j in 1:w){
#    minc <- sum(pvalues[j,Kc[,cb]] <=min(ct)) + 1
#    maxc <- sum(pvalues[j,Kc[,cb]] <=max(ct))
#    
#    if(is.null(shift)){shift = 0}
#    lambdaA[j] <- min(sapply(c(minc:maxc), function(x) quant[j, Kc[   sort.list(pvalues[j,Kc[,cb]])[x] ,cb   ] ]/x  ))
#    #lamb <- min(lamb,  betaquantsS[, combs2[   sort.list(pvmatr.uns[j,combs2[,c]]),c]][j,a]/a  )
    
    
#  }
  
  
  # nRej_Perm <- rejPerm(pvalues,min(ct))
  # quantRej_min <- sort(nRej_Perm)[ceiling((1-alpha)*w)]
  # if(quantRej_min>0){
  #   lambdaE <- min(lambdaE, (min(ct) + shift)*m/(quantRej_min*alpha))}
  
#  return(lambdaA)
#}

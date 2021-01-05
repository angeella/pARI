###########################Lambda Optimal############################

#For every permutation, I compute the best possible lambda. Then I take the alpha w +1quantile due to discreteness.

#The idea with the min and max cut-offs is that you only take into account the p-values between 
#these cut-offs. So the smaller the interval (1utoff, mutoff) is, the more power you will have. 
#The interval (1utoff, mutoff) corresponds with the set \Gamma in the paper.

#pvalues = ordered pvalues 
#ct= vector of cutoffs
#family = Beta or Simes

lambdaOptR <- function(pvalues, family, alpha, delta){
  #pvalues matrix with dimensions variables times permutations
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  library(Rcpp)
  source("src/rowSortC.cpp")
  if(is.unsorted(pvalues[1,])){pvalues = rowSortC(pvalues)}
  #implement set of threshold
  
  l <- c()
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  if(is.null(delta) ){delta = 0}
  for(j in 1:w){
    if(family =="simes"){
      
      lambda <- ((m-delta)*(pvalues[j,1:m]))/((c(1:m)-delta)*alpha)
    }
    if(family == "beta"){
      lambda <- pbeta(q =pvalues[j,c(1:m)],shape1 = c(1:m),shape2 =m+1-c(1:m))
      
    }
    
    if(family == "finner"){
      lambda <- (pvalues[j,c(1:m)]*(m - c(1:m) + delta) )/ (alpha * (c(1:m) - delta - (pvalues[j,c(1:m)]* (c(1:m) - 1))))
      
    }
    
    
    #   if(family =="higher.criticism"){
    
    #     lambda <- (sqrt(m)*((c(1:m)/m) - pvalues[j,c(1:m)]))/(sqrt(pvalues[j,c(1:m)]*(1-pvalues[j,c(1:m)])))
    
    #HigherCriticism <- function(lambda,pvalue){
    
    # m <- length(pvalue)
    #  i <- c(1:m)
    #  (2*i + lambda^2 - sqrt((2*i + lambda^2)^2 - 4*i^2 * (m + lambda^2)/m))/(2*(m+ lambda^2)) - pvalue
    
    #}
    #lambda <- sapply(c(1:m), function(x) rootSolve::uniroot.all(HigherCriticism,pvalue = pvalues[j,x],lower=0, upper=1000))
    #lambda <- unlist(lambda)
    #}
    l[j] <- min(lambda)
    
    
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
#    1 <- sum(pvalues[j,Kc[,cb]] <=min(ct)) + 1
#    m <- sum(pvalues[j,Kc[,cb]] <=max(ct))
#    
#    if(is.null(shift)){shift = 0}
#    lambdaA[j] <- min(sapply(c(1:m), function(x) quant[j, Kc[   sort.list(pvalues[j,Kc[,cb]])[x] ,cb   ] ]/x  ))
#    #lamb <- min(lamb,  betaquantsS[, combs2[   sort.list(pvmatr.uns[j,combs2[,c]]),c]][j,a]/a  )


#  }


# nRej_Perm <- rejPerm(pvalues,min(ct))
# quantRej_min <- sort(nRej_Perm)[ceiling((1-alpha)*w)]
# if(quantRej_min>0){
#   lambdaE <- min(lambdaE, (min(ct) + shift)*m/(quantRej_min*alpha))}

#  return(lambdaA)
#}

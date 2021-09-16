
#' @importFrom Rmpfr mpfr  
#' @importFrom Rmpfr .mpfr_erange_set
#' @importFrom Rmpfr .mpfr_erange
###########################Lambda Optimal############################

#For every permutation, I compute the best possible lambda. Then I take the alpha w +1quantile due to discreteness.

#The idea with the min and max cut-offs is that you only take into account the p-values between 
#these cut-offs. So the smaller the interval (1utoff, mutoff) is, the more power you will have. 
#The interval (1utoff, mutoff) corresponds with the set \Gamma in the paper.

#pvalues = ordered pvalues 
#ct= vector of cutoffs
#family = Beta or Simes

lambdaOptR <- function(pvalues, family, alpha, delta, m = NULL){
  #pvalues matrix with dimensions variables times permutations
  family_set <- c("simes", "aorc", "beta", "higher.criticism")
  #sourceCpp("src/rowSortC.cpp")
  
  family <- match.arg(tolower(family), family_set)
  if(is.unsorted(pvalues[1,])){pvalues = rowSortC(pvalues)}
  #implement set of threshold
  
  l <- c()
  w <- dim(pvalues)[1]
  if(is.null(m)){m <- dim(pvalues)[1]}
  if(is.null(delta) ){delta = 0}
  for(j in 1:w){
    if(family =="simes"){
      
      lambda <- ((m-delta)*(pvalues[j,1:m]))/((c(1:m)-delta))
    }
    if(family == "beta"){
     # lambda <- pbeta(q =pvalues[j,c(1:m)],shape1 = c(1:m),shape2 =m+1-c(1:m))
      md <- c(1:m)/(m+1)
      sdd <- sqrt((c(1:m)*(m+1-c(1:m)))/((m+1)^2 * (m+2)))
      .mpfr_erange_set(value = (1-2^-52)*.mpfr_erange(c("min.emin","max.emax")))
      lambda <- pnorm(q = mpfr(pvalues[j, c(1:m)], 128), mean = md, 
                              sd = sdd)
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
  
  lambdaE <- stats::quantile(, alpha, type = 1)
  

  
  return(lambdaE)
}




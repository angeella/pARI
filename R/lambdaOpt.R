#' @title Lambda calibration
#' @description \code{lambdaOpt} provides the optimal lambda calibration parameter used in the critical vector.
#' @usage lambdaOpt(pvalues,family,alpha, delta, step.down, max.step)
#' @param pvalues pvalues matrix with dimensions equal to the number of variables times the number of permutations.
#' @param family by default \code{family="simes"}. Choose a family of confidence envelopes to compute the critical vector from \code{"simes"}, \code{"AORC"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param alpha alpha level.
#' @param delta by default \code{delta = 0}. Do you want to consider sets with at least delta size?
#' @param step.down by default \code{step.down = FALSE}. If you want to compute the lambda calibration parameter using the step down approach put TRUE.
#' @param max.step by default \code{max.step = 10}. Maximum number of steps for the step down approach
#' @author Angela Andreella
#' @return lambda parameter
#' @export
#' @importFrom stats pbeta

lambdaOpt <- function(pvalues, family, alpha, delta, step.down = FALSE, max.step = 10){
  #pvalues matrix with dimensions variables times permutations
  family_set <- c("simes", "aorc", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  lambdaE <- lambdaCalibrate(X = pvalues, alpha = alpha, delta = delta, family = family)
  #lambdaE <- lambdaOpt1(pvalues= t(pvalues), alpha = alpha, delta = delta, family = family)
  
  if(step.down){
    convergence <- FALSE
    cv0 <- criticalVector(pvalues = pvalues, family = family, alpha = alpha, lambda = lambdaE, delta = delta)
    no_rej <- which(pvalues[,1] > min(cv0))
    it <- 1
    lambda <- c()
    while(!convergence && it < max.step){
      lambda[it] <- lambdaCalibrate(X = pvalues[no_rej,], alpha = alpha, delta = delta, family = family)
      cv1 <- criticalVector(pvalues = pvalues[no_rej,], family = family, alpha = alpha, lambda = lambda[it], delta = delta)
      no_rej_new <- which(pvalues[,1] > min(cv1))
      
      if(all(no_rej_new %in% no_rej)){
        convergence <- TRUE
        lambdaE <- max(lambdaE,lambda[it])
      }else{
        no_rej <- no_rej_new
      }
    }
  }
  
 # if(is.unsorted(pvalues[,1])){pvalues = rowSortC(pvalues)}
  #implement set of threshold
  
  
  
#  l <- c()
#  w <- dim(pvalues)[1]
#  m <- dim(pvalues)[2]
#  if(is.null(delta) ){delta = 0}
#  for(j in 1:w){
#    minc <- sum(pvalues[j,] <=min(ct)) + 1
#    maxc <- sum(pvalues[j,] <=max(ct))
#    if(family =="simes"){
      
#      minc = minc + delta
 #     lambda <- ((m-delta)*(pvalues[j,minc:maxc]))/((c(minc:maxc)-delta)*alpha)
      #lambda <- (m*(pvalues[j,minc:maxc] + shift))/(c(minc:maxc))
#    }
#    if(family == "beta"){
#      lambda <- pbeta(q =pvalues[j,c(minc:maxc)],shape1 = c(minc:maxc),shape2 =m+1-c(minc:maxc))
      
#    }
    
#    if(family == "finner"){
#      minc = minc + delta
      #lambda <- (pvalues[j,c(minc:maxc)]*(m - c(minc:maxc) + delta) )/ (alpha * (c(minc:maxc) - delta - (pvalues[j,c(minc:maxc)]* (c(minc:maxc) - 1))))
 #     lambda <- (pvalues[j,c(minc:maxc)]*(m - 1) )/ (alpha * (c(minc:maxc) - delta) * (1 - pvalues[j,c(minc:maxc)]))
      
 #   }
    
    
 #   if(family =="higher.criticism"){
      
 #     lambda <- (sqrt(m)*((c(minc:maxc)/m) - pvalues[j,c(minc:maxc)]))/(sqrt(pvalues[j,c(minc:maxc)]*(1-pvalues[j,c(minc:maxc)])))
      
      #HigherCriticism <- function(lambda,pvalue){
        
      # m <- length(pvalue)
      #  i <- c(1:m)
      #  (2*i + lambda^2 - sqrt((2*i + lambda^2)^2 - 4*i^2 * (m + lambda^2)/m))/(2*(m+ lambda^2)) - pvalue
        
      #}
      #lambda <- sapply(c(minc:maxc), function(x) rootSolve::uniroot.all(HigherCriticism,pvalue = pvalues[j,x],lower=0, upper=1000))
      #lambda <- unlist(lambda)
#}
#    l[j] <- min(lambda)
    
    
#  }
  
#  lambdaE <- sort(l)[floor(alpha*w)+1]

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

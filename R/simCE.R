#' @title Simulation Confidence Envelope
#' @description Simulation Confidence Envelope using the gaussianSamples function from the package sansSouci
#' @usage simCE(m,n,pi0,B,rho,SNR,alpha,delta,ct,family)
#' @param m Number of variables, i.e. hypotheses
#' @param n Number of observations
#' @param pi0 Proportion of true null hypotheses
#' @param B number of permutations
#' @param rho Level of equi-correlation between pairs of variables
#' @param SNR Signal to noise ratio
#' @param alpha significance level
#' @param delta which size set do you want to consider at least?
#' @param ct set of threshold
#' @param family family of confidence bound
#' @return simulation CE
#' @export
#' @importFrom sansSouci gaussianSamples

simCE <- function(m,n,pi0,B,rho,SNR,alpha,delta,ct,family){
  
  sim <- gaussianSamples(m, rho, n, pi0, SNR = SNR, prob = 1)
  X <- sim$X #rows represent variables
  out <- signTest(X,B)
  P <- t(cbind(out$pv, out$pv_H0)) #rows represent the permutations  
  P_ord <- rowSortC(P)
  p <- P[1,]
  
  lambda <- lambdaOpt(P_ord, family = family, ct = ct, alpha = alpha, delta = delta)
  cv <- cv(pvalues=P_ord, family= family, alpha = alpha, delta = delta, lambda = lambda, ct = ct)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  #h <- hI(p, cv)
  
  #cvh <- sapply(c(1:length(p)), function(x) ((x * alpha * lambda)/h)- shift)
  
  out <- as.data.frame(cbind(cv, p, t(P), 
                             rep(n = length(cv), lambda),
                             rep(n = length(cv), B),
                             rep(n = length(cv), pi0),
                             rep(n = length(cv), SNR),
                             rep(n = length(cv), rho),
                             rep(n = length(cv), delta)))
  colnames(out) <- c("CriticalValue", 
                     "Pvalues", 
                     sapply(c(1:(B+1)),function(x) paste0("B_unsort",x)), 
                     "lambda",
                     "B", "pi0", "SNR", "rho", "delta"
  )
  
  return(out = out)
  
}

simCEVectorize <- Vectorize(simCE, vectorize.args = c("pi0", "SNR", "rho"),SIMPLIFY = FALSE)

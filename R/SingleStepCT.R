#' @title Single Step Permutation-based method
#' @description Performs single step method based on confidence envelope constructed by sign flipping 
#' @usage SingleStepCT(X,ct, ix, alpha, family, delta= NULL, B)
#' @param X data matrix (rows represent variables, columns observations), the first rows are the raw pvalues.
#' @param ct set of thresholds of interest
#' @param alpha alpha level
#' @param family family of confidence envelopes considered in order to find the critical values. Now, we implement the Simes ones and the Beta
#' @param delta do you want to consider at least delta size set?
#' @param ix set of hypothesis of interest
#' @param B number of permutations
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{discoveries} number of discoveries in the set selected, \code{pvalues} raw pvalues
#' @export

SingleStepCT <- function(X,ct, ix, alpha, family, delta= NULL, B){

  out <- signTest(X, B = B)
  P <- t(cbind(out$pv, out$pv_H0))
  P_ord <- rowSortC(P)
  p <- P[1,]

  lambda <- lambdaOpt(P_ord, family = family, ct = ct, alpha = alpha, delta = delta)
  cv <- cv(pvalues=P_ord, family= family, alpha = alpha, delta = delta, lambda = lambda, ct = ct)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  #h <- hI(p, cv)
  
  #cvh <- sapply(c(1:length(p)), function(x) ((x * alpha * lambda)/h)- shift)
  
  discoveries <- dI(ix,cv,p)
  TDP <- discoveries/length(ix)
  
  return(list(discoveries = discoveries, TDP = TDP, pvalues = p, cv = cv))
  
}

cvhPerm <- function(praw, alpha, shift, family, lambda){
  #p <- pvalues[1,]

  cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  h <- hI(praw, cv)
  
  cvh <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/h)- shift)
  
  return(cvh)
  
}

discoveriesPerm <- function(out, ix){
  praw <- out[[3]]
  cv <- out[[4]]
  discoveries <- dI(ix,cv,praw)
  TDP <- discoveries/length(ix)
  
  return(list(discoveries = discoveries, TDP = TDP, praw = praw))
  
}

#a<-SingleStepCT(pvalues=pvalues,ct=c(0.001,0.01), ix = c(1:4000), alpha = 0.1, shift = 0, family='Simes', lambda =1)
#a[[1]]
#discoveries(hommel(p = pvalues[1,],simes = TRUE),alpha=0.1,ix = c(1:4000))
#I have two discoveries more.


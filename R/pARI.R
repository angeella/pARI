#' @title Single Step permutation-based method
#' @description Performs single step method based on confidence envelope constructed by sign flipping 
#' @usage pARI(X,ix, alpha, alternative, family, delta, B, pvalues)
#' @param X data matrix rows represent variables, columns observations, the first rows are the raw pvalues.
#' @param ix set of hypothesis of interest
#' @param alpha alpha level, default 0.1
#' @param alternative character referring to the alternative hypothesis, \code{"two.sided"}, \code{"greater"} or \code{"lower"}. Default is \code{"two.sided"}
#' @param family family of confidence envelopes considered in order to find the critical values. Now, we implement the Simes ones and the Beta
#' @param delta do you want to consider at least delta size set?
#' @param B number of permutations
#' @param pvalues matrix pvalues instead of data with dimensions hypotheses times permutations, default NULL
#' @author Angela Andreella
#' @return Returns a list with the following objects discoveries number of discoveries in the set selected, pvalues raw pvalues
#' @export

pARI <- function(X= NULL, ix, alpha = 0.1, alternative = "two.sided", family = "simes", delta= 0, B = 1000, pvalues = NULL){
  if(is.null(pvalues)){
    out <- signTest(X, B = B, alternative = alternative)
    P <- cbind(out$pv, out$pv_H0)
  }else{
    P <- pvalues
  }

  p <- P[,1]

  lambda <- lambdaOpt(P, family = family, alpha = alpha, delta = delta)
  cv <- cv(pvalues=P, family= family, alpha = alpha, delta = delta, lambda = lambda)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  #h <- hI(p, cv)
  
  #cvh <- sapply(c(1:length(p)), function(x) ((x * alpha * lambda)/h)- shift)
  
  discoveries <- dI(ix,cv,p)
  TDP <- discoveries/length(ix)
  
  return(list(discoveries = discoveries, TDP = TDP, pvalues = p, cv = cv))
  
}

#cvhPerm <- function(praw, alpha, shift, family, lambda){
  #p <- pvalues[1,]

#  cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
  
  #Compute the largest size of a set of hyp not rejected by our local test
#  h <- hI(praw, cv)
  
#  cvh <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/h)- shift)
  
#  return(cvh)
  
#}




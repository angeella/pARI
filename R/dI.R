#' @title Lower bound true discoveries
#' @description Calculates (1-alpha) lower confidence bounds for the set-wise of false null hypotheses
#' @usage dI(ix,cv,praw)
#' @param ix set-wise of hypotheses considered
#' @param cv vector of critical values
#' @param praw raw pvalues
#' @param iterative if \code{iterative = TRUE}, the iterative iterative method for improvement of confidence envelopes is applied. 
#' @param approx if \code{iterative = TRUE} and you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb if \code{approx = TRUE}, you must decide how many large random subcollection (level of approximation) considered.
#' @author Angela Andreella
#' @return Returns the lower confidence bound
#' @export

dI <- function(ix,cv,praw, iterative, approx, ncomb){
  
#  u <- sapply(c(1:length(ix)), function(x) 1 - x + sum(praw[ix] <= cv[x]))
#  d <- max(u)
  if(iterative){
  d <- permDiscoveriesIt(ix, cv, praw, approx, ncomb)
  }else{
  d <- permDiscoveries(ix,cv,praw)
  }
  return(d)
}


#######################Largest size of hyp not rejected by the local test###############

#Lemma 8 admissible paper GS 2019

#pvalues = vector of raw pvalues
#cv = critical values of our local test.

#hI <- function(praw, cv){
  
#  praw_sort <- sort(praw)
#  cv_ord <- sort(cv)
#  m <- length(praw_sort)
#  size <- c()
#  n <- function(N){
#    cond <- sum(sapply(c(1:N), function(x) praw_sort[m - N + x] > cv_ord[x]))
#  }
#  nV <- Vectorize(n, vectorize.args = "N")  
#  h <- nV(c(1:m))
  
#  h <- max(h[h==c(1:m)])
#  return(h)
#}

#ix = set of hyp
#cv = set of critical vector computed in h
#praw = pvalues raw
#h = from hI







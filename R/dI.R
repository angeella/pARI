#' @title Lower bound true discoveries
#' @description Calculates (1-alpha) lower confidence bounds for the set-wise of false null hypotheses
#' @usage dI(ix, cv, pvalues, iterative, approx, ncomb, family, alpha, delta)
#' @param ix set-wise of hypotheses considered
#' @param cv vector of critical values
#' @param pvalues pvalues matrix with dimensions variables times permutations
#' @param iterative if \code{iterative = TRUE}, the iterative iterative method for improvement of confidence envelopes is applied. 
#' @param approx if \code{iterative = TRUE} and you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb if \code{approx = TRUE}, you must decide how many large random subcollection (level of approximation) considered.
#' @param family if \code{iterative = TRUE} specify the family of confidence envelopes.
#' @param alpha if \code{iterative = TRUE} specify the alpha level.
#' @param delta if \code{iterative = TRUE} specify the delta level.
#' @export
#' @author Angela Andreella
#' @return Returns the lower confidence bound

dI <- function(ix, cv, pvalues, iterative, approx, ncomb, family, alpha, delta){
#  u <- sapply(c(1:length(ix)), function(x) 1 - x + sum(praw[ix] <= cv[x]))
#  d <- max(u)

  d <- permDiscoveries(ix = ix, cv = cv, praw = pvalues[,1])
  print(d)
  if(iterative){
    d_seq <- c()
    d_seq[1] <- d
    print(d_seq)
    it <- 1
    dist <- Inf
    while(dist !=0) {
      if(approx == TRUE){
        Kcomb <- replicate(ncomb, sample(ix,size =length(ix) - d_seq[it], replace=FALSE), 
                           simplify="matrix")
      }else{
        Kcomb <- combn(ix, length(ix) - d_seq[it]) 
      }
      #Create complementry set: combinations + all not in ix
      m <- dim(pvalues)[1]
      R <- which(!(c(1:m) %in% ix))
      print("2")
      lambda_kc <- sapply(c(1:ncomb), function(x) {
        if(is.matrix(Kcomb)){
          Kc <- Kcomb[,x]
        }else{Kc <- Kcomb[x]}
        Kc <- unique(c(Kc, R))
        P_Kc <- matrix(pvalues[Kc,], nrow= length(Kc), ncol = dim(pvalues)[2])
        lambdaCalibrate(X = P_Kc, family = family, alpha = alpha, delta = delta)
      })
      print("3")
      print(lambda_kc)
      lambda <- max(lambda_kc,lambda)
      print(lambda)
      cv <- criticalVector(pvalues= pvalues, family= family, 
                           alpha = alpha, delta = delta, lambda = lambda)
      
      d_seq[it +1] <- permDiscoveries(ix = ix,cv = cv, praw = pvalues[,1])
      
      dist <- d_seq[it] - d_seq[it+1]
      it <- it + 1
      
    }
    print("4")
    print(d_seq)
    # B_est <- min(Bt)
    
    d <- max(d_seq)
    
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







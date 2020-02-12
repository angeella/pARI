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

dI <- function(ix,cv,praw){
  
  u <- sapply(c(1:length(ix)), function(x) 1 - x + sum(praw[ix] <= cv[x]))
  d <- max(u)
  return(d)
}





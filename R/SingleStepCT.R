###############################################Single Step Method#################################

#pvalues = matrix of pvalues, where the first rows are the raw pvalues. We use testByRandomization by sansSouci
#package. It computes the one sample t test using sign-flipping
#ct = set of thresholds of interest
#alpha
#ix = set of hypothesis of interest
#shift = shift the lines of critical values in order to consider larger pvalues
#family = family of confidence envelopes considered in order to find the critical values. Now, we implement the 
#Simes ones and the Beta


SingleStepCT <- function(pvalues,ct, ix, alpha, shift, family, lambda, ...){

  p <- pvalues[1,]
  if(is.null(shift)){shift = 0}
  cv <- cv(pvalues=P_ord, family= family, alpha = alpha, shift = 7, lambda = lambda, ct = ct)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  #h <- hI(p, cv)
  
  #cvh <- sapply(c(1:length(p)), function(x) ((x * alpha * lambda)/h)- shift)
  
  discoveries <- dI(ix,cv,p)
  
  return(list(discoveries, p))
  
}

cvhPerm <- function(praw, alpha, shift, family, lambda){
  #p <- pvalues[1,]

  cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  h <- hI(praw, cv)
  
  cvh <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/h)- shift)
  
  return(cvh)
  
}

discoveriesPerm <- function(praw, ix, cvh){
  
  discoveries <- dI(ix,cvh,praw)
  
  return(list(discoveries, praw))
  
}

#a<-SingleStepCT(pvalues=pvalues,ct=c(0.001,0.01), ix = c(1:4000), alpha = 0.1, shift = 0, family='Simes', lambda =1)
#a[[1]]
#discoveries(hommel(p = pvalues[1,],simes = TRUE),alpha=0.1,ix = c(1:4000))
#I have two discoveries more.


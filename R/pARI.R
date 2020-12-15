#' @title Single Step permutation-based method
#' @description Performs single step method based on confidence envelope constructed by sign flipping 
#' @usage pARI(X,ix, alpha, family, delta, B, pvalues, test.type, complete, ...)
#' @param X data matrix rows represent variables, columns observations, the first rows are the raw pvalues.
#' @param ix set of hypothesis of interest. It can be a vector having same lenght of the variables or the indices if one set is considered.
#' @param alpha alpha level, default 0.1
#' @param family family of confidence envelopes considered in order to find the critical values. Now, we implement the Simes ones and the Beta
#' @param delta do you want to consider at least delta size set?
#' @param B number of permutations
#' @param pvalues matrix pvalues instead of data with dimensions hypotheses times permutations, default NULL
#' @param test.type "one_sample", i.e., one sample, or "two_samples", i.e., two samples  t-tests?
#' @param complete default FALSE. If TRUE the sets of critical vectors and raw pvalues are returned.
#' @param ... Futher arguments, see details.
#' @seealso \code{\link{signTest}} \code{\link{permTest}}
#' @author Angela Andreella
#' @return Returns a list with the following objects discoveries number of discoveries in the set selected, pvalues raw pvalues
#' @export

pARI <- function(X= NULL, ix, alpha = 0.05, family = "simes", delta = 0, B = 1000, pvalues = NULL, test.type = "one_sample", complete = FALSE, ...){
 
  #Add different design then two and one
  #Check for error
   if(is.null(pvalues)){
     if(test.type == "one_sample"){
       out <- signTest(X, B = B, ...)
     }else{
       out <- permTest(X, B = B, ...)
     }
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
  
  if(length(ix) == nrow(X)){
  levels_ix <- unique(ix)  
  discoveries <- c()
  ixX <- list()
  TDP <- c()
    for(i in 1:length(levels_ix)){
      ixX[[i]] <- which(ix == levels_ix[i])
      discoveries[i] <- dI(ixX[[i]],cv,p)
      TDP[i] <- discoveries[i]/length(ixX[[i]])
      if(!is.null(rownames(X))){
      ixX[[i]] <- rownames(X)[ixX[[i]]]
        }
    }

  }else{
    discoveries <- dI(ix,cv,p)
    TDP <- discoveries/length(ix)
    if(!is.null(rownames(X))){
      ixX <- rownames(X)[ix]
    }else{
      ixX <- ix
    }
  }
  

  if(complete){
    return(list(discoveries = discoveries, TDP = TDP, ix = ixX, pvalues = p, cv = cv))
  }else{
    return(list(discoveries = discoveries, TDP = TDP, ix = ixX))
    
  }
  
}

#cvhPerm <- function(praw, alpha, shift, family, lambda){
  #p <- pvalues[1,]

#  cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
  
  #Compute the largest size of a set of hyp not rejected by our local test
#  h <- hI(praw, cv)
  
#  cvh <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/h)- shift)
  
#  return(cvh)
  
#}




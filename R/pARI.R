#' @title Permutation-based All-Resolutions Inference
#' @description The main function for single step All-Resolutions Inference (ARI) method based on critical vectors constructed by permutations. 
#' @usage pARI(X,ix, alpha, family, delta, B, pvalues, test.type, complete, ...)
#' @param X data matrix where rows represent variables, and columns the observations.
#' @param ix set of hypotheses of interest. It can be a vector having the same length as the variables or the indices of the variables of interest if only one set is considered.
#' @param alpha alpha level, default 0.05.
#' @param family by default \code{family="simes"}. Choose a family of confidence envelopes to compute the critical vector from \code{"simes"}, \code{"finner"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param delta by default \code{delta = 0}. Do you want to consider sets with at least delta size?
#' @param B by default \code{B = 1000}. Number of permutations.
#' @param pvalues by default \code{pvalues = NULL}.  Matrix of pvalues used instead of the data matrix having dimensions equal to the number of hypotheses times the number of permutations.
#' @param test.type by default \code{test.type = "one_sample"}. Choose a type of tests among \code{"one_sample"}, i.e., one-sample t-test, or \code{"two_samples"}, i.e., two-samples t-tests.
#' @param complete by default \code{complete = FALSE}. If \code{TRUE} the sets of critical vectors and the raw pvalues are returned.
#' @param ... Futher parameters.
#' @seealso The type of tests implemented: \code{\link{signTest}} \code{\link{permTest}}.
#' @author Angela Andreella
#' @return by default returns a list with the following objects: \code{discoveries}: lower bound for the number of true discoveries in the set selected, \code{ix}: selected variables. If \code{complete = TRUE} the raw \code{pvalues} and \code{cv} critical vector are returned.
#' @export
#' @references For the general framework of All-Resolutions Inference see:
#' 
#' Goeman, Jelle J., and Aldo Solari. "Multiple testing for exploratory research." Statistical Science 26.4 (2011): 584-597.
#'
#' For permutation-based All-Resolutions Inference see:
#' 
#' Andreella, Angela, et al. "Permutation-based true discovery proportions for fMRI cluster analysis." arXiv preprint arXiv:2012.00368 (2020).
#' 
#' @examples
#' datas <- simulateData(pi0 = 0.8, m = 1000, n = 30, power = 0.9, rho = 0,set.seed = 123)
#' out <- pARI(X = datas, ix = c(1:200),test.type = "one_sample")
#' out

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
  cv <- criticalVector(pvalues=P, family= family, alpha = alpha, delta = delta, lambda = lambda)
  
  #Compute the largest size of a set of hyp not rejected by our local test
  #h <- hI(p, cv)
  
  #cvh <- sapply(c(1:length(p)), function(x) ((x * alpha * lambda)/h)- shift)
  
  if(length(ix) == nrow(P)){
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
    discoveries <- dI(ix = ix,cv = cv,praw = p)
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




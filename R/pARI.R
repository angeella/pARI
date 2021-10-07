#' @title Permutation-based All-Resolutions Inference
#' @description The main function for All-Resolutions Inference (ARI) method based on critical vectors constructed by permutations. 
#' @usage pARI(X= NULL, ix, alpha = 0.05, family = "simes", delta = 0, B = 1000, pvalues = NULL, 
#' test.type = "one_sample", complete = FALSE, clusters = FALSE, iterative = FALSE, 
#' approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...)
#' @param X data matrix where rows represent the \code{m} variables and columns the \code{n} observations.
#' @param ix set of hypotheses of interest. It can be a vector having the same length as the variables indicating the respective cluster, in this case you must put \code{clusters = TRUE}, or the indices of the variables of interest if only one set/cluster is considered.
#' @param alpha alpha level.
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.#' @param alpha alpha level.
#' @param delta delta value. Do you want to consider sets with at least delta size? By default \code{delta = 0}. 
#' @param B Number of permutations, by default \code{B = 1000}. 
#' @param pvalues by default \code{pvalues = NULL}.  Matrix of pvalues used instead of the data matrix having dimensions equal to the number of hypotheses times the number of permutations.
#' @param test.type by default \code{test.type = "one_sample"}. Choose a type of tests among \code{"one_sample"}, i.e., one-sample t-test, or \code{"two_samples"}, i.e., two-samples t-tests.
#' @param complete by default \code{complete = FALSE}. If \code{TRUE} the sets of critical vectors and the raw pvalues are returned.
#' @param clusters if \code{ix} indicates the clusters/sets must be \code{TRUE}
#' @param iterative if \code{iterative = TRUE}, the iterative iterative method for improvement of confidence envelopes is applied. Default is \code{FALSE}.
#' @param approx if \code{iterative = TRUE} and you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb if \code{approx = TRUE}, you must decide how many random subcollections (level of approximation) considered.
#' @param step.down by default \code{step.down = FALSE}. If you want to compute the lambda calibration parameter using the step-down approach put \code{TRUE}.
#' @param max.step by default \code{max.step = 10}. Maximum number of steps for the step down approach, so useful when \code{step.down = TRUE}.
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

pARI <- function(X= NULL, ix, alpha = 0.05, family = "simes", delta = 0, B = 1000, pvalues = NULL, test.type = "one_sample", complete = FALSE, clusters = FALSE, 
                 iterative = FALSE, approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10,...){
 
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

  lambda <- lambdaOpt(P, family = family, alpha = alpha, delta = delta, step.down = step.down, max.step = max.step)
  cv <- criticalVector(pvalues=P, family= family, alpha = alpha, delta = delta, lambda = lambda)
  
  if(length(ix) == nrow(P) & clusters == TRUE){
  levels_ix <- unique(ix)  
  discoveries <- c()
  ixX <- list()
  TDP <- c()
    for(i in 1:length(levels_ix)){
      ixX[[i]] <- which(ix == levels_ix[i])
      discoveries[i] <- dI(ixX[[i]],cv,P, iterative, approx, ncomb, family, alpha, delta)
      TDP[i] <- discoveries[i]/length(ixX[[i]])
      if(!is.null(rownames(X))){
      ixX[[i]] <- rownames(X)[ixX[[i]]]
        }
    }

  }else{
    discoveries <- dI(ix = ix,cv = cv,pvalues = P, iterative, approx, ncomb, family, alpha, delta)
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





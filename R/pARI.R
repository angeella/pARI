#' @title Permutation-based All-Resolutions Inference
#' @description The main function for All-Resolutions Inference (ARI) method based on critical vectors constructed 
#' using the p-values permutation distribution. The function computes simultaneous lower bounds for the number of true discoveries 
#' for each set of hypotheses specified in \code{ix} controlling family-wise error rate.
#' @usage pARI(X= NULL, ix, alpha = 0.05, family = "simes", delta = 0, B = 1000, pvalues = NULL, 
#' test.type = "one_sample", complete = FALSE, clusters = FALSE, iterative = FALSE, 
#' approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...)
#' @param X data matrix where rows represent the \eqn{m} variables and columns the \eqn{n} observations.
#' @param ix numeric vector which expresses the set of hypotheses of interest. It can be a vector with length equals \eqn{m} indicating the corresponding cluster for each variable,
#' (in this case, you must put \code{clusters = TRUE}), or a vector containing the position indices of the variables of interest if only one set/cluster of hypotheses is considered.
#' @param alpha numeric value in `[0,1]`. It expresses the alpha level to control the family-wise error rate.
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.#' @param alpha alpha level.
#' @param delta numeric value. It expresses the delta value, please see the references. Default to 0. 
#' @param B numeric value. Number of permutations, default to 1000. 
#' @param pvalues matrix of pvalues with dimensions \eqn{m \times B} used instead of the data matrix \code{X}. Default to @NULL.
#' @param test.type character string. Choose a type of tests among \code{"one_sample"}, i.e., one-sample t-test, or \code{"two_samples"}, i.e., two-samples t-tests. Default \code{"one_sample"}.
#' @param complete Boolean value. If \code{TRUE} the sets of critical vectors and the raw pvalues are returned. Default @FALSE. 
#' @param clusters Boolean value. If \code{ix} indicates many clusters/sets must be \code{TRUE}. Default @FALSE.
#' @param iterative Boolean value. If \code{iterative = TRUE}, the iterative method for improvement of confidence envelopes is applied. Default @FALSE.
#' @param approx Boolean value. Default @TRUE. If you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb Numeric value. If \code{approx = TRUE}, you must decide how many random subcollections (level of approximation) considered.
#' @param step.down Boolean value. Default @FALSE If you want to compute the lambda calibration parameter using the step-down approach put \code{TRUE}.
#' @param max.step Numeric value. Default to 10. Maximum number of steps for the step down approach, so useful when \code{step.down = TRUE}.
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
#' datas <- simulateData(pi0 = 0.8, m = 1000, n = 30, power = 0.9, rho = 0,seed = 123)
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





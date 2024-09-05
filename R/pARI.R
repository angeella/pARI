#' @title Permutation-based All-Resolutions Inference
#' @description The main function for All-Resolutions Inference (ARI) method based on the critical vector constructed 
#' using the \eqn{p}-values permutation distribution. The function computes simultaneous lower bounds for the number of true discoveries 
#' for each set of hypotheses specified in \code{ix} controlling family-wise error rate at level \code{alpha}.
#' @usage pARI(X= NULL, ix, alpha = 0.05, family = "simes", delta = 0, B = 1000, pvalues = NULL, 
#' test.type = "one_sample", complete = FALSE, clusters = FALSE, iterative = FALSE, 
#' approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...)
#' @param X Data matrix where rows represent the \eqn{m} variables and columns the \eqn{n} observations.
#' @param ix Numeric vector which expresses the set of hypotheses of interest. It can be a vector with length equals \eqn{m} indicating the corresponding cluster for each variable,
#' (in this case, you must put \code{clusters = TRUE}), or a vector containing the position indices of the variables of interest if only one set/cluster of hypotheses is considered.
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @param family String character. Name of the family confidence envelope to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"}, \code{"higher.criticism"}, and \code{"power"}.
#' Default to "simes".
#' @param delta Numeric value. \eqn{\delta} value. Please see the reference below. Default to 0. 
#' @param B Numeric value. Number of permutations, default to 1000. 
#' @param pvalues Matrix of \eqn{p}-values with dimensions \eqn{m \times B} where \eqn{m} is the number of variables 
#' and \eqn{B} the number of permutations used instead of the data matrix \code{X}. Default to NULL.
#' @param test.type Character string. Choose a type of tests among \code{"one_sample"}, i.e., one-sample t-tests, or \code{"two_samples"}, i.e., two-samples t-tests. Default \code{"one_sample"}.
#' @param complete Boolean value. If \code{TRUE} the sets of critical vectors and the raw \eqn{p}-values are returned. Default to \code{FALSE}. 
#' @param clusters Boolean value. If \code{ix} indicates many clusters/sets must be \code{TRUE}. Default @FALSE.
#' @param iterative Boolean value. If \code{iterative = TRUE}, the iterative method is applied (computationally demanding). Default to \code{FALSE}. Please see the reference below.
#' @param approx Boolean value. Default to \code{TRUE}. If you are analyzing high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time. Please see the reference below.
#' @param ncomb Numeric value. If \code{approx = TRUE}, you must decide how many random sub collections (level of approximation) considered. Default to 100.
#' @param step.down Boolean value. Default to \code{FALSE} If you want to compute the lambda calibration parameter using the step-down approach put \code{TRUE}. Please see the reference below.
#' @param max.step Numeric value. Default to 10. Maximum number of steps for the step down approach, so useful when \code{step.down = TRUE}.
#' @param ... Further arguments
#' @seealso The type of tests implemented: \code{\link{signTest}} \code{\link{permTest}}.
#' @author Angela Andreella
#' @return by default returns a list with the following objects: 
#' \describe{
#'  \item{discoveries}{lower bound for the number of true discoveries in the set selected}
#'  \item{ix}{selected variables} 
#' }
#' If \code{complete = TRUE} the raw \code{pvalues} and \code{cv} critical vector are also returned.
#' @export
#' @references For the general framework of All-Resolutions Inference see:
#' 
#' Goeman, Jelle J., and Aldo Solari. "Multiple testing for exploratory research." Statistical Science 26.4 (2011): 584-597.
#'
#' For permutation-based All-Resolutions Inference see:
#' 
#' Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.
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
  P <- ifelse(is.na(P),1, P)
  p <- P[,1]

  lambda <- lambdaOpt(pvalues = P, family = family, alpha = alpha, delta = delta, step.down = step.down, max.step = max.step)
  cv <- criticalVector(pvalues=P, family= family, alpha = alpha, delta = delta, lambda = lambda)
  
  if(length(ix) == nrow(P) & clusters == TRUE){
  levels_ix <- unique(ix)  
  discoveries <- c()
  ixX <- list()
  TDP <- c()
    for(i in 1:length(levels_ix)){
      ixX[[i]] <- which(ix == levels_ix[i])
      discoveries[i] <- dI(ix = ixX[[i]],cv = cv,pvalues = P, 
                           iterative = iterative, approx = approx, ncomb = ncomb, 
                           family = family, alpha = alpha, delta = delta)
      TDP[i] <- discoveries[i]/length(ixX[[i]])
      if(!is.null(rownames(X))){
      ixX[[i]] <- rownames(X)[ixX[[i]]]
        }
    }

  }else{
    discoveries <- dI(ix = ix,cv = cv,pvalues = P, 
                      iterative = iterative, approx = approx, ncomb = ncomb, 
                      family = family, alpha = alpha, delta = delta)
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





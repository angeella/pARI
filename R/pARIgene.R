#' @title Permutation-based All-Resolutions Inference for Gene Expression Data
#' @description This function computes the lower bound for the  number of true 
#' discoveries within each cluster (pathways) of Gene Expression Data.
#' @usage pARIgene(X= NULL, pathways, alpha = 0.05, family = "simes", delta = 0, 
#' B = 1000, test.type = "one_sample", complete = FALSE, iterative = FALSE, 
#' approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...)
#' @param X data matrix where rows represent the \eqn{m} variables and columns 
#' the \eqn{n} observations.
#' @param pathways list of pathways where names indicates the name of the pathway.
#' @param alpha numeric value in `[0,1]`. It expresses the alpha level to control 
#' the family-wise error rate.
#' @param family string character. Choose a family of confidence envelopes to 
#' compute the critical vector from \code{"simes"}, \code{"aorc"}, \code{"beta"} 
#' and \code{"higher.criticism"}.#' @param alpha alpha level.
#' @param delta numeric value. It expresses the delta value, 
#' please see the references. Default to 0. 
#' @param B numeric value. Number of permutations, default to 1000. 
#' @param test.type character string. Choose a type of tests among 
#' \code{"one_sample"}, i.e., one-sample t-test, or \code{"two_samples"}, 
#' i.e., two-samples t-tests. Default \code{"one_sample"}.
#' @param complete Boolean value. If \code{TRUE} the sets of critical vectors 
#' and the raw pvalues are returned. Default @FALSE. 
#' @param iterative Boolean value. If \code{iterative = TRUE}, the iterative 
#' method for improvement of confidence envelopes is applied. Default @FALSE.
#' @param approx Boolean value. Default @TRUE. If you are treating high 
#' dimensional data, we suggest to put \code{approx = TRUE} to speed up 
#' the computation time.
#' @param ncomb Numeric value. If \code{approx = TRUE}, you must decide 
#' how many random subcollections (level of approximation) considered.
#' @param step.down Boolean value. Default @FALSE If you want to compute 
#' the lambda calibration parameter using the step-down approach put \code{TRUE}.
#' @param max.step Numeric value. Default to 10. Maximum number of steps for 
#' the step down approach, so useful when \code{step.down = TRUE}.
#' @param ... Futher parameters.
#' @seealso The type of tests implemented: \code{\link{signTest}} 
#' \code{\link{permTest}}.
#' @author Angela Andreella
#' @return by default returns a list with the following objects: 
#' \code{discoveries}: lower bound for the number of true discoveries 
#' in the set selected, \code{ix}: selected variables. I
#' f \code{complete = TRUE} the raw \code{pvalues} and \code{cv} 
#' critical vector are returned.
#' @export
#' @references For the general framework of All-Resolutions Inference see:
#' 
#' Goeman, Jelle J., and Aldo Solari. "Multiple testing for exploratory research.
#' " Statistical Science 26.4 (2011): 584-597.
#'
#' For permutation-based All-Resolutions Inference see:
#' 
#' Andreella, Angela, et al. "Permutation-based true discovery proportions for 
#' fMRI cluster analysis." arXiv preprint arXiv:2012.00368 (2020).
#' 

pARIgene <- function(X= NULL, pathways, alpha = 0.05, family = "simes", 
                     delta = 0, B = 1000, test.type = "one_sample", 
                     complete = FALSE, iterative = FALSE, approx = TRUE, 
                     ncomb = 100, step.down = FALSE, max.step = 10,...){
  
  #Add different design then two and one
  #Check for error
  if(test.type == "one_sample"){
      out <- signTest(X, B = B, ...)
    }else{
      out <- permTest(X, B = B, ...)
    }
    P <- cbind(out$pv, out$pv_H0)

  P <- ifelse(is.na(P),1, P)
  p <- P[,1]
  
  lambda <- lambdaOpt(pvalues = P, family = family, alpha = alpha, 
                      delta = delta, step.down = step.down, max.step = max.step)
  cv <- criticalVector(pvalues=P, family= family, alpha = alpha, delta = delta, 
                       lambda = lambda)
  
  discoveries <- c()
  ixX <- list()
  TDP <- c()
  for(i in seq(pathways)){
    ixX[[i]] <- which(rownames(X) %in% pathways[[i]])
    discoveries[i] <- dI(ix = ixX[[i]],cv = cv,pvalues = P, 
                           iterative = iterative, approx = approx, ncomb = ncomb, 
                           family = family, alpha = alpha, delta = delta)
    TDP[i] <- discoveries[i]/length(ixX[[i]])
    if(!is.null(rownames(X))){
      ixX[[i]] <- rownames(X)[ixX[[i]]]
    }
  }
    
  full <- dI(ix = seq(nrow(P)),cv = cv,pvalues = P, 
             iterative = iterative, approx = approx, ncomb = ncomb, 
             family = family, alpha = alpha, delta = delta)
  

  out <- data.frame(discoveries = c(discoveries, full), 
                    TDP = c(round(TDP,5), round(full/nrow(P), 5)), 
                    pathways = c(names(pathways), "full_genes"))
  
  if(complete){
    return(list(out = out,
                pvalues = p, cv = cv))
  }else{
    return(out)
    
  }
  
}
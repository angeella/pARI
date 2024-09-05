#' @title Permutation-based All-Resolutions Inference for Gene Expression Data
#' @description This function computes the lower bound for the  number of true 
#' discoveries within each cluster (pathways) of Gene Expression Data.
#' @usage pARIgene(X= NULL, pathways, alpha = 0.05, family = "simes", delta = 0, 
#' B = 1000, test.type = "one_sample", complete = FALSE, iterative = FALSE, 
#' approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...)
#' @param X Data matrix where rows represent the \eqn{m} variables and columns the \eqn{n} observations.
#' @param pathways List of pathways where names indicates the name of the pathway.
#' @param alpha Numeric value in `[0,1]`. \eqn{\alpha} level to control the family-wise error rate. Default to 0.05.
#' @param family String character. Name of the family confidence envelope to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"}, \code{"higher.criticism"}, and \code{"power"}.
#' Default to "simes".
#' @param delta Numeric value. \eqn{\delta} value. Please see the reference below. Default to 0. 
#' @param B Numeric value. Number of permutations, default to 1000. 
#' @param test.type Character string. Choose a type of tests among \code{"one_sample"}, i.e., one-sample t-tests, or \code{"two_samples"}, i.e., two-samples t-tests. Default \code{"one_sample"}.
#' @param complete Boolean value. If \code{TRUE} the sets of critical vectors and the raw \eqn{p}-values are returned. Default to \code{FALSE}. 
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
#' Goeman, Jelle J., and Aldo Solari. "Multiple testing for exploratory research.
#' " Statistical Science 26.4 (2011): 584-597.
#'
#' For permutation-based All-Resolutions Inference see:
#' 
#' Andreella, A., Hemerik, J., Finos, L., Weeda, W., & Goeman, J. (2023). Permutation-based true discovery proportions for functional magnetic resonance imaging cluster analysis. Statistics in Medicine, 42(14), 2311-2340.


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
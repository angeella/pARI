#' @title Permutatation-based one sample t-test
#' @description Performs sign-flipped one-sample t-tests.
#' @usage signTest(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F)
#' @param X data matrix where rows represent the \code{m} variables and columns the \code{n} observations.
#' @param B numeric value, number of permutations to be performed, including the identity. Default is 1000.
#' @param alternative character string referring to the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}). 
#' @param seed numeric value, specify seed. Default is 1234.
#' @param mask 3D array of logicals (i.e. \code{TRUE/FALSE} in/out of the brain). Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that none of the voxels have to be excluded.
#' @param rand logical value, \code{TRUE} to compute p-values by permutation distribution.
#' @author Angela Andreella
#' @return Returns a list with the following objects:
#' \describe{ 
#'   \item{Test}{Vector of \eqn{m} observed one-sample t-tests}
#'   \item{Test_H0}{Matrix with dimensions \eqn{m \times B} of permuted one-sample t-tests}
#'   \item{pv}{Vector of \eqn{m} observed p-values} 
#'   \item{pv_H0}{Matrix with dimensions \eqn{m \times B} of permuted p-values}}
#' @export
#' @importFrom stats pnorm
#' @importFrom RNifti readNifti
#' @importFrom matrixStats rowRanks
#' @examples 
#' X <- matrix(rnorm(100*20), ncol=20)
#' out <- oneSample(X = X, alternative = "two.sided")

signTest <- function(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F){
  
  alternative_set <- c("two.sided", "greater", "lower")
  
  if(!is.null(seed)){set.seed(seed)}else{set.seed(1234)}
  if(!is.null(mask)){
    if(is.character(mask)){mask = readNifti(mask)}
    X <- X[which(mask==1),]
    
  }
  
  alternative <- match.arg(tolower(alternative), alternative_set)
  
  ## number of obeservation
  n <- ncol(X)
  # number of variables
  m <- nrow(X)
  
  ## Observed test statistics
  rowV <- rowVariance(X)
  Test <- ifelse(rowV==0,0, rowMeans(X)/(sqrt((rowV)/n)))
  
  ## Test statistics under H0
  
  Test_H0 <- signFlip(X,B-1)
  
  Test_H0 <- ifelse(is.na(Test_H0), 0 , Test_H0)
  
  if(!rand){
    pv <- switch(alternative, 
                 "two.sided" = 2*(pt(abs(Test), df = n-1, lower.tail=FALSE)),
                 "greater" = pt(Test, df = n-1, lower.tail=FALSE),
                 "lower" = 1-pt(Test, df = n-1, lower.tail=FALSE))
    
    pv_H0 <- switch(alternative, 
                    "two.sided" = 2*(pt(abs(Test_H0), df = n-1,  lower.tail=FALSE)),
                    "greater" = pt(Test_H0, df = n-1,  lower.tail=FALSE),
                    "lower" = 1-pt(Test_H0, df = n-1,  lower.tail=FALSE)) 
  }else{

    Test_matrix <- cbind(Test, Test_H0)
    pv_matrix <- switch(alternative, 
                 "two.sided" = rowRanks(-abs(Test_matrix)) / ncol(Test_matrix),
                 "greater" = rowRanks(-Test_matrix) / ncol(Test_matrix),
                 "lower" = rowRanks(Test_matrix) / ncol(Test_matrix))
    
    pv <- pv_matrix[, 1]
    pv_H0 <- pv_matrix[, 2:(B)]
  }
  
  out <- list(Test = Test, Test_H0 = Test_H0, pv = pv, pv_H0 = pv_H0)

  return(out)
}
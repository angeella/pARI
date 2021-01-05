#' @title Sign-Flipping Test
#' @description Performs sign-flipping, i.e. permutation, one-sample t-tests
#' @usage signTest(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F)
#' @param X data where rows represents the variables and columns the observations
#' @param B by default \code{B = 1000}. Number of permutations.
#' @param alternative a character string referring to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param seed by default \code{seed=1234}. Integer value specifying the seed.
#' @param mask 3D array of locicals (i.e. \code{TRUE/FALSE} in/out of the brain). Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @param rand by default \code{rand = FALSE}. 
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Test} observed one sample t-test, \code{Test_H0} Test statistics under H0, \code{pv} observed p-values, \code{pv_H0} p-values under H0.
#' @export
#' @importFrom stats pnorm
#' @importFrom RNifti readNifti
#' @importFrom matrixStats rowRanks


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
                 "greater" = pt(Test,  lower.tail=FALSE),
                 "lower" = 1-pt(Test,  lower.tail=FALSE))
    
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
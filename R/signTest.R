#' @title Sign-Flipping Test
#' @description Performs sign-flipping, i.e. permutation, one-sample t-tests
#' @usage signTest(X, B = 1000, alternative = "two.sided", seed = NULL)
#' @param X data where rows represents the variables and columns the observations
#' @param B number of permutations to perform, default is 1000.
#' @param alternative character referring to the alternative hypothesis, "two.sided", "greater" or "lower". Default is "two.sided"
#' @param seed specify seed, default is 1234.
#' @param mask mask to apply, it can be nii format or path
#' @param rand logical. Should p values computed by permutation distribution?
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Test} observed one sample t-test, \code{Test_H0} Test statistics under H0, \code{pv} observed p-values, \code{pv_H0} p-values under H0
#' @export
#' @importFrom stats pnorm
#' @importFrom RNifti readNifti
#' @importFrom matrixStats rowRanks

alternative_set <- c("two.sided", "greater", "lower")

signTest <- function(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F){
  

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
  #rowV <- ifelse(rowV==0,.Machine$double.xmin, rowV)
  #Test <- rowMeans(X)/(sqrt((rowV)/n))
  Test <- ifelse(rowV==0,0, rowMeans(X)/(sqrt((rowV)/n)))
  
  ## Test statistics under H0
  
  Test_H0 <- signFlip(X,B)
  #T0_m <- meanBySignFlipping(X,B)
  #T0_v <- varBySignFlipping(X,B)
  #T0_v <- ifelse(T0_v==0,.Machine$double.xmin, T0_v)
  #Test_H0 <- T0_m/ sqrt((T0_v)/n)
  Test_H0 <- ifelse(is.na(Test_H0), 0 , Test_H0)
  
  if(!rand){
    pv <- switch(alternative, 
                 "two.sided" = 2*(pnorm(abs(Test), lower.tail=FALSE)),
                 "greater" = pnorm(Test, lower.tail=FALSE),
                 "lower" = 1-pnorm(Test, lower.tail=FALSE))
    
    pv_H0 <- switch(alternative, 
                    "two.sided" = 2*(pnorm(abs(Test_H0), lower.tail=FALSE)),
                    "greater" = pnorm(Test_H0, lower.tail=FALSE),
                    "lower" = 1-pnorm(Test_H0, lower.tail=FALSE))
  }else{

    Test_matrix <- cbind(Test, Test_H0)
    pv_matrix <- switch(alternative, 
                 "two.sided" = rowRanks(-abs(Test_matrix)) / (B+1),
                 "greater" = rowRanks(-Test_matrix) / (B+1),
                 "lower" = rowRanks(Test_matrix) / (B+1))
    
    pv <- pv_matrix[, 1]
    pv_H0 <- pv_matrix[, 2:(B+1)]
  }
  
  out <- list(Test = Test, Test_H0 = Test_H0, pv = pv, pv_H0 = pv_H0)

  return(out)
}
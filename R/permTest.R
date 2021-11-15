#' @title Permutation Test
#' @description Performs permutation-based two-sample t-tests.
#' @usage permTest(X, B = 1000, alternative = "two.sided", seed = NULL, 
#' mask = NULL, rand = FALSE, label = NULL)
#' @param X data matrix where rows represent the \eqn{m} variables and columns the \eqn{n} observations.
#' @param B numeric value. Number of permutations, default to 1000. 
#' @param alternative character string. It refers to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param seed integer value. If you want to specify the seed. Default to @NULL
#' @param mask NIfTI file or character string. 3D array of logical values (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @param rand Boolean value. Default @FALSE. If \code{rand = TRUE}, the p-values are computed by \code{\link{rowRanks}}.
#' @param label numeric/character vector. Labels of the observations, if \code{NULL} the columns's name are considered. Default @NULL. 
#' @author Angela Andreella
#' @return Returns a list with the following objects: 
#' - \code{Test}: vector with length equals \eqn{m}. Observed two-samples t-tests, one for each \eqn{m} variable, 
#' - \code{Test_H0}: matrix with dimensions \eqn{m \times B-1}. Test statistics under H0,
#' - \code{pv}: vector with length equals \eqn{m}. observed p-values, one for each \eqn{m} variable,
#' - \code{pv_H0} matrix with dimensions \eqn{m \times B-1}. P-values under H0.
#' @export
#' @importFrom stats pnorm
#' @importFrom RNifti readNifti
#' @importFrom matrixStats rowRanks
#' @importFrom stats pt


permTest <- function(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = FALSE, label = NULL){
  
  alternative_set <- c("two.sided", "greater", "lower")
  if(!is.null(seed)){set.seed(seed)}else{set.seed(1234)}
  if(!is.null(mask)){
    if(is.character(mask)){mask = readNifti(mask)}
    X <- X[which(mask==1),]
    
  }
  
  alternative <- match.arg(tolower(alternative), alternative_set)
  
  if(is.null(label)){label <- colnames(X)}
  
  #id <- unique(label)
  label <- factor(label)
  levels(label) <- c(0,1)
  ## number of obeservation
  n <- ncol(X)
  # number of variables
  m <- nrow(X)
  
  ## Observed test statistics
  id <- levels(label)
  n1 <- sum(label==id[1])
  n2 <- sum(label==id[2])
  rowV1 <- rowVariance(X[,label == id[1]])
  rowV2 <-rowVariance(X[,label == id[2]])
  rowM1 <- rowMeans(X[,label == id[1]])
  rowM2 <-rowMeans(X[,label == id[2]])
  #pooled.var <- ((n1 - 1)* rowV1 + (n2 - 1)* rowV2)/ (n1 + n2 - 2)
  pooled.var <- (rowV1/n1 + rowV2/n2)
  #Test <- (rowM1 - rowM2)/sqrt(pooled.var)
  Test <- (rowM1 - rowM2)/sqrt(pooled.var)
  Test <- ifelse(is.na(Test), 0 , Test)
  ## Test statistics under H0

  Test_H0 <- permT(as.matrix(X),B-1,label)
  Test_H0 <- ifelse(is.na(Test_H0), 0 , Test_H0)
  
  if(!rand){
    
    gdl <- ((rowV1/n1 + rowV2/n2)^2)/((((rowV1/n1)^2)/(n1-1))+(((rowV2/n2)^2)/(n2-1)))
  #  gdl <- n1 + n2 - 2
    pv <- switch(alternative, 
                 "two.sided" = 2*(pt(abs(Test), df = gdl, lower.tail=FALSE)),
                 "greater" = pt(Test, df = gdl, lower.tail=FALSE),
                 "lower" = 1-pt(Test, df = gdl, lower.tail=FALSE))
    
    pv_H0 <- switch(alternative, 
                    "two.sided" = 2*(pt(abs(Test_H0), df = gdl, lower.tail=FALSE)),
                    "greater" = pt(Test_H0, df = gdl, lower.tail=FALSE),
                    "lower" = 1-pt(Test_H0, df = gdl, lower.tail=FALSE))
  }else{
    
    Test_matrix <- cbind(Test, Test_H0)
    pv_matrix <- switch(alternative, 
                        "two.sided" = rowRanks(-abs(Test_matrix)) / (B),
                        "greater" = rowRanks(-Test_matrix) / (B),
                        "lower" = rowRanks(Test_matrix) / (B))
    
    pv <- pv_matrix[, 1]
    pv_H0 <- pv_matrix[, 2:(B)]
  }
  
  out <- list(Test = Test, Test_H0 = Test_H0, pv = pv, pv_H0 = pv_H0)
  
  return(out)
}
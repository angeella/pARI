#' @title Permutation Test
#' @description Performs permutation tests, i.e. permutation, two-sample t-tests
#' @usage permTest(X, B, alternative, seed, mask, rand, label)
#' @param X data where rows represents the variables and columns the observations, the columns's name define the categories
#' @param B number of permutations to perform, default is 1000.
#' @param alternative character referring to the alternative hypothesis, "two.sided", "greater" or "lower". Default is "two.sided"
#' @param seed specify seed, default is 1234.
#' @param mask mask to apply, it can be nii format or path
#' @param rand logical. Should p values computed by permutation distribution?
#' @param label labels of the observations, if NULL the columns's name are considered, default NULL
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Test} observed one sample t-test, \code{Test_H0} Test statistics under H0, \code{pv} observed p-values, \code{pv_H0} p-values under H0
#' @export
#' @importFrom stats pnorm
#' @importFrom RNifti readNifti
#' @importFrom matrixStats rowRanks
#' @importFrom stats pt


permTest <- function(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F, label = NULL){
  
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
  id <- c(0,1)
  n1 <- sum(label==id[1])
  n2 <- sum(label==id[2])
  rowV1 <- rowVariance(X[,label == id[1]])
  rowV2 <-rowVariance(X[,label == id[2]])
  rowM1 <- rowMeans(X[,label == id[1]])
  rowM2 <-rowMeans(X[,label == id[2]])
  #pooled.var <- ((n1 - 1)* rowV1 + (n2 - 1)* rowV2)/ (n1 + n2 - 2)
  pooled.var <- (rowV1/n1 + rowV2/n2)
  #Test <- (rowM1 - rowM2)/sqrt(pooled.var)
  Test <- (rowM1 - rowM2)/sqrt(pooled.var * (1/n1 + 1/n2))
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
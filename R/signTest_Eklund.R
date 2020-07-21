
alternative_set <- c("two.sided", "greater", "lower")

signTest_Eklund <- function(X, B = 1000, alternative = "two.sided", seed = NULL, mask = NULL, rand = F){
  
  library(matrixStats)
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
  
  scoresFlip <- X
  for(i in 1:ncol(scoresFlip)){
    
    scoresFlip[,i] <- sapply(c(1:nrow(scoresFlip)), function(x) scoresFlip[x,i] %*% (rbinom(1, 1, 0.5)*2 - 1))
    
  }
  rowV <- rowVars(scoresFlip) + rowMeans(scoresFlip)^2
  Test <- ifelse(is.na(rowV) | rowV == 0,0, rowMeans(scoresFlip)/(sqrt((rowV))))
  ## Test statistics under H0
  
  Test_H0 <- signFlip_Eklund(scoresFlip,B-1)
  #T0_m <- meanBySignFlipping(X,B)
  #T0_v <- varBySignFlipping(X,B)
  #T0_v <- ifelse(T0_v==0,.Machine$double.xmin, T0_v)
  #Test_H0 <- T0_m/ sqrt((T0_v)/n)
  Test_H0 <- ifelse(is.na(Test_H0), 0 , Test_H0)
  
  if(!rand){
    pv <- switch(alternative, 
                 "two.sided" = 2*(pt(abs(Test), df = n-1, lower.tail=FALSE)),
                 "greater" = pt(Test, df = n-1, lower.tail=FALSE),
                 "lower" = 1-pt(Test, df = n-1, lower.tail=FALSE))
    
    pv_H0 <- switch(alternative, 
                    "two.sided" = 2*(pt(abs(Test_H0), df = n-1, lower.tail=FALSE)),
                    "greater" = pt(Test_H0, df = n-1, lower.tail=FALSE),
                    "lower" = 1-pt(Test_H0, df = n-1, lower.tail=FALSE))
  }else{
    
    Test_matrix <- cbind(Test, Test_H0)
    pv_matrix <- switch(alternative, 
                        "two.sided" = rowRanks(-abs(Test_matrix)) / (B+1),
                        "greater" = rowRanks(-Test_matrix) / (B+1),
                        "lower" = rowRanks(Test_matrix) / (B+1))
    
    pv <- pv_matrix[, 1]
    pv_H0 <- pv_matrix[, 2:(B)]
  }
  
  out <- list(Test = Test, Test_H0 = Test_H0, pv = pv, pv_H0 = pv_H0)
  
  return(out)
}
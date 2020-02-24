oneSample<- function(X,alternative = c("two.sided", "less", "greater")){
  
  rowVars <- function (x,na.rm = TRUE) 
  {
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
  }
  n <- ncol(X)
  m <- nrow(X)
  alternative <- match.arg(alternative)
  
  T = rowMeans(X)/(sqrt((rowVars(X))/n))
  p <- switch(alternative, 
              #"two.sided" = 2*(1 - pnorm(abs(T))),
              "two.sided" = 2*(pnorm(abs(T), lower.tail=FALSE)),
              #"greater" = 1 - pnorm(T),
              "greater" = pnorm(T, lower.tail=FALSE),
              #"less" = pnorm(T))
              "less" = 1-pnorm(T, lower.tail=FALSE))
  
  res <- list(T = T, p = p)
  
  return(res)
}
  

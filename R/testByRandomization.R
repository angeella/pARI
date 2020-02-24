testByRandomization <- function(X, B, 
                                alternative = c("two.sided", "less", "greater"),
                                rowTestFUN = rowWelchTests,
                                rand.p.value = FALSE, seed = NULL){
  rowVars <- function (x,na.rm = TRUE) 
  {
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
  }
  
  alternative <- match.arg(alternative)
  ## sanity checks
  n <- ncol(X)
  m <- nrow(X)

  ## observed test statistics and p-values
  #T <- rowSums(X)/sqrt(2*n)
  T <- rowMeans(X)/(sqrt((rowVars(X))/n))
  p <- switch(alternative, 
              #"two.sided" = 2*(1 - pnorm(abs(T))),
              "two.sided" = 2*(pnorm(abs(T), lower.tail=FALSE)),
              #"greater" = 1 - pnorm(T),
              "greater" = pnorm(T, lower.tail=FALSE),
              #"less" = pnorm(T))
              "less" = 1-pnorm(T, lower.tail=FALSE))
    ## test statistics under H0
    #T0 <- signFlip(X,B)
    T0_m <- meanBySignFlipping(X,B)
    T0_v <- varBySignFlipping(X,B)
    T0 <- T0_m/ sqrt((T0_v)/n)
    p0 <- switch(alternative, 
                 #"two.sided" = 2*(1 - pnorm(abs(T0))),
                 "two.sided" = 2*(pnorm(abs(T0), lower.tail=FALSE)),
                 #"greater" = 1 - pnorm(T0),
                 "greater" = pnorm(T0, lower.tail=FALSE),
                 "less" = 1-pnorm(T0, lower.tail=FALSE))
    #"less" = pnorm(T0))
    res <- list(T = T, T0 = T0, p = p, p0 = p0)

  return(res)
}
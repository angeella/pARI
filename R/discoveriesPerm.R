#' @title Lower bound for the number of true discoveries
#' @description Internal function.
#' @usage discoveriesPerm(out, ix) 
#' @param out output from the \code{\link{pARI}} function.
#' @param ix numeric vector. It refers to the set-wise hypotheses considered. 
#' @author Angela Andreella
#' @return discoveriesPerm


discoveriesPerm <- function(out, ix){
  praw <- out[[3]]
  cv <- out[[4]]
  discoveries <- dI(ix = ix,cv = cv,praw = praw)
  TDP <- discoveries/length(ix)
  
  return(list(discoveries = discoveries, TDP = TDP, praw = praw))
  
}


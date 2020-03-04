#' @title discoveriesPerm
#' @description compute discoveriesPerm
#' @usage discoveriesPerm(out, ix) 
#' @param out out
#' @param ix ix
#' @author Angela Andreella
#' @return discoveriesPerm
#' @export

discoveriesPerm <- function(out, ix){
  praw <- out[[3]]
  cv <- out[[4]]
  discoveries <- dI(ix,cv,praw)
  TDP <- discoveries/length(ix)
  
  return(list(discoveries = discoveries, TDP = TDP, praw = praw))
  
}

#a<-SingleStepCT(pvalues=pvalues,ct=c(0.001,0.01), ix = c(1:4000), alpha = 0.1, shift = 0, family='Simes', lambda =1)
#a[[1]]
#discoveries(hommel(p = pvalues[1,],simes = TRUE),alpha=0.1,ix = c(1:4000))
#I have two discoveries more.
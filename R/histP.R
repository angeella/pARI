#' @title Histogram Pvalues
#' @description create histogram p-values
#' @usage histP(copes, alternative, mask, zstat, path, name)
#' @param copes image copes instead of pvalues, default NULL
#' @param alternative alternative
#' @param mask mask
#' @param zstat zstat
#' @param path path
#' @param name name
#' @param method method to compute the zstat
#' @return Returns hist
#' @export
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics hist
#' @importFrom grDevices rgb
#' @importFrom graphics legend

histP <- function(copes =NULL, alternative, mask, zstat = NULL, path = getwd(), name = "hist",method=NULL){
  

  if(is.null(copes) & is.null(zstat)){stop('Please insert observed pvalues or copes images')}
  
  if(!is.null(copes)){
    
    if(is.null(mask)){stop('please insert the mask as character path or Nifti image')}
    if(is.character(mask)){mask = readNifti(mask)}
    if(!is.list(copes)){stop("Please insert the list of copes as list class object")}
    
    img_dims <- c(91,  109 , 91)
    img <- array(NA, c(img_dims, length(copes)))
    
    for (sid in 1:length(copes)) {  
      img[,,,sid] <- copes[[sid]]
      
    }
    
    scores <- matrix(img,nrow=(91*109*91),ncol=length(copes))
    scores <- scores[which(mask==1),]
    res <- signTest(X=scores, B = 1000,alternative = alternative) #variables times number of permutation
    
    pvalues <- cbind(res$pv,res$pv_H0)
    pvalues = t(pvalues)
    rm(res)
    rm(scores)
    rm(copes)
    rm(img)
    
    #pvalues_ord <- rowSortC(pvalues)
    
    
  }
    if(!is.null(zstat) & is.character(zstat)){
      zstat = readNifti(zstat)
      zstat = get_array(zstat)
      zstat = zstat[which(mask==1)]
      praw <- switch(alternative, 
                   "two.sided" = 2*(pnorm(abs(zstat), lower.tail=FALSE)),
                   "greater" = pnorm(zstat, lower.tail=FALSE),
                   "less" = 1-pnorm(zstat, lower.tail=FALSE))
    }else{
      zstat = get_array(zstat)
      zstat = zstat[which(mask==1)]
      praw <- switch(alternative, 
                     "two.sided" = 2*(pnorm(abs(zstat), lower.tail=FALSE)),
                     "greater" = pnorm(zstat, lower.tail=FALSE),
                     "less" = 1-pnorm(zstat, lower.tail=FALSE))
    }
  
    png(paste0(path,"/", name, ".png")) 
    hist(pvalues[1,], col = rgb(0, 1, 1, 0.5), xlab = "p-values", main = "P-values Histogram")
    for(i in 2:nrow(pvalues)){
      
      hist(pvalues[i,],col=rgb(0, 0, 0, 0.15),add=TRUE)
      
    }
    hist(pvalues[1,], col = rgb(0, 1, 1, 0.5),add=TRUE)
    if(!is.null(praw)){
      hist(praw, col = rgb(1, 0, 0, 0.5),add=TRUE)
      legend('top',c("Observed", "Null", paste0(method)), col= c(rgb(0, 1, 1, 0.5),rgb(0, 0, 0, 0.15),rgb(1, 0, 0, 0.5)),lwd =4)
    }else{
      legend('top',c("Observed", "Null"), col= c(rgb(0, 1, 0, 0.5),rgb(0, 0, 0, 0.15)),lwd =4)
      
    }
    dev.off()
  
  
  
}

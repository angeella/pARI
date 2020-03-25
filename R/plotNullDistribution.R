#' @title Plot Permutation Pvalues
#' @description create plot permutation pvalues with corresponding critical vectors
#' @usage plotNullDistribution(P,family,alpha, ct, path, name, delta,copes,mask)
#' @param P permutation matrix pvalues where rows correspond to the permutation, default NULL
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is NULL
#' @param alpha type I error allowed to construct the critical family, default 0.1
#' @param ct set of threshold to construct the critical family, default c(0,1)
#' @param path selected the path where the plot will be saved, if NULL the current working directory is used
#' @param name name plot pdf, if NULL plot is used
#' @param delta for the family critical values, default NULL
#' @param copes image copes instead of pvalues, default NULL
#' @param mask mask
#' @return Returns plot null distribution with critical value curve and observed pvalues in red
#' @export
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @importFrom graphics legend

plotNullDistribution <- function(P=NULL,family="simes",alpha = 0.1, ct = c(0,1), path = getwd(), name = "plot", delta = NULL,copes=NULL,mask=NULL){
  
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  fam_match <- function(x) {match.arg(tolower(x), family_set)}
  if(!is.null(family)){family <- unlist(lapply(family, fam_match))}
  if(is.null(copes) & is.null(P)){stop('Please insert pvalues matrix or copes images')}
  
  if(!is.null(P) & is.unsorted(P[1,])){pvalues_ord <- rowSortC(P)}
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
    res <- signTest(X=scores, B = 1000,alternative = "two.sided") #variables times number of permutation
    
    pvalues <- cbind(res$pv,res$pv_H0)
    pvalues = t(pvalues)
    rm(res)
    rm(scores)
    rm(copes)
    rm(img)
    
    pvalues_ord <- rowSortC(pvalues)
    
    
  }

  if(is.null(family)){
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[1,], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:nrow(pvalues_ord)){
      
      lines(pvalues_ord[i,],col='black',type="l")
      
    }
    lines(pvalues_ord[1,], lwd =2, col= 'red')
    dev.off()
  }else{
    lcv <- function(family,delta=NULL, cols = "blue"){
      lambdaO <- lambdaOpt(pvalues = pvalues_ord,family=family,ct=ct,alpha=alpha, delta = delta)
      cvO<- cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda = lambdaO, delta = delta)
      lines(cvO, lwd =2, col= cols)
    }
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    lcvV <- Vectorize(lcv,vectorize.args = c("family", "delta", "cols"))
    cols = rainbow(length(family))
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[1,], type = 'l', col = ' green', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:nrow(pvalues_ord)){
      
      lines(pvalues_ord[i,],col='black',type="l")
      
    }
    lines(pvalues_ord[1,], lwd =2, col= 'green')
    #lines(cvO, col= 'blue', lwd =2)
    mapply(lcv, family, delta, cols)
    family <- firstup(family)
    legend('top',legend=c(sapply(c(1:length(family)), 
                                 function(x) as.expression(bquote(~ .(family[x]) ~ delta == .(delta[x]) ))), 
                          " Observed Pvalues"), col= c(cols, "green"),lwd =2)
    
    dev.off()
  }
  
  
}

#' @title Plot Permutation Pvalues
#' @description create plot permutation pvalues with corresponding critical vectors
#' @usage plotNullDistribution(P,family,alpha, path, name, delta,copes,mask, alternative, rand, B)
#' @param P permutation matrix pvalues where rows correspond to the variables, default NULL
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is NULL. Vectors accepted
#' @param alpha type I error allowed to construct the critical family, default 0.1
#' @param path selected the path where the plot will be saved, if NULL the current working directory is used
#' @param name name plot pdf, if NULL plot is used
#' @param delta for the family critical values, default 0. Vectors accepted
#' @param copes image copes instead of pvalues, default NULL
#' @param mask mask
#' @param alternative character referring to the alternative hypothesis, "two.sided", "greater" or "lower". Default is "two.sided"
#' @param rand logical. Should p values computed by permutation distribution?
#' @param B number of permutations to perform, default is 1000, if matrix p-value is computed directly.
#' @return Returns plot null distribution with critical value curve and observed pvalues in red
#' @export
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @importFrom graphics legend

plotNullDistribution <- function(P=NULL,family="simes",alpha = 0.1, path = getwd(), name = "plot", delta = 0,copes=NULL,mask=NULL, alternative = "two.sided", rand = F, B = 1000){
  
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  fam_match <- function(x) {match.arg(tolower(x), family_set)}
  alternative_set <- c("two.sided", "greater", "lower")
  alternative <- match.arg(tolower(alternative), alternative_set)
  if(!is.null(family)){family <- unlist(lapply(family, fam_match))}
  if(is.null(copes) & is.null(P)){stop('Please insert pvalues matrix or copes images')}
  
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
    res <- signTest(X=scores, B = B,alternative = alternative, rand = rand) #variables times number of permutation
    
    P <- cbind(res$pv,res$pv_H0)
    rm(res)
    rm(scores)
    rm(copes)
    rm(img)
  
    
  }

  if(is.null(family)){
    if(is.unsorted(P[,1])){pvalues_ord <- colSortC(P)}else{pvalues_ord <- P}
    
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[,1], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:ncol(pvalues_ord)){
      
      lines(pvalues_ord[,i],col='black',type="l")
      
    }
    lines(pvalues_ord[,1], lwd =2, col= 'red')
    dev.off()
  }else{
    if(is.unsorted(P[,1])){pvalues_ord <- colSortC(P)}else{pvalues_ord <- P}
    
    lcv <- function(family,delta=NULL, cols = "blue"){
      lambdaO <- lambdaOpt(pvalues = P,family=family,alpha=alpha, delta = delta)
      cvO<- cv(pvalues = P, family = family, alpha = alpha, lambda = lambdaO, delta = delta)
      lines(cvO, lwd =2, col= cols)
    }
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    lcvV <- Vectorize(lcv,vectorize.args = c("family", "delta", "cols"))
    cols = rainbow(length(family))
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[,1], type = 'l', col = ' red', lty = "dashed", xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:ncol(pvalues_ord)){
      
      lines(pvalues_ord[,i],col='black',type="l")
      
    }
    lines(pvalues_ord[,1], lwd =2, col= 'red', lty = "dashed")
    #lines(cvO, col= 'blue', lwd =2)
    mapply(lcv, family, delta, cols)
    family <- firstup(family)
    legend('top',legend=c(sapply(c(1:length(family)), 
                                 function(x) as.expression(bquote(~ .(family[x]) ~ delta == .(delta[x]) ))), 
                          " Observed Pvalues"), col= c(cols, "red"),lwd =2, lty =c(rep("solid", length(family)), "dashed"))
    
    dev.off()
  }
  
  
}

#Plot Permutation Pvalues
##
# @param P permutation matrix pvalues where rows correspond to the permutation
# @param family Do you want to plot some cv family? It can be one or more
# @param alpha type I error allowed to construct the critical family
# @param ct set of threshold to construct the critical family
# @param path selected the path where the plot will be saved, if NULL the current working directory is used
# @param name name plot pdf, if NULL plot is used
# @param delta for the family critical values

plotNullDistribution <- function(P,family=NULL,alpha = 0.1, ct = c(0,1), path = getwd(), name = "plot", delta = NULL){
  
  pvalues_ord <- rowSortC(P)

  if(is.null(family)){
    pdf(paste0(path,"/", name, ".pdf")) 
    plot(pvalues_ord[1,], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:nrow(pvalues_ord)){
      
      lines(pvalues_ord[i,],col='black',type="l")
      
    }
    lines(pvalues_ord[1,], lwd =2, col= 'red')
    dev.off()
  }else{
    lambdaO <- lambdaOpt(pvalues = pvalues_ord,family=family,ct=ct,alpha=alpha, delta = delta)
    cvO <- cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda = lambdaO, delta = delta)

    pdf(paste0(path,"/", name, ".pdf")) 
    plot(pvalues_ord[1,], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:nrow(pvalues_ord)){
      
      lines(pvalues_ord[i,],col='black',type="l")
      
    }
    lines(pvalues_ord[1,], lwd =2, col= 'red')
    lines(cvO, col= 'blue', lwd =2)
    dev.off()
  }
  
  
}
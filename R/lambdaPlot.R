#' @title Lambda plot
#' @description compute lambda calibration plot
#' @usage lambdaPlot(copes,family,ct,alpha, delta,P, mask, path, name)
#' @param copes list of copes, default NULL
#' @param family family
#' @param alpha alpha
#' @param ct set threshold
#' @param delta delta
#' @param P matrix of pvalues where rows indicate the permutations, default NULL
#' @param mask mask
#' @param path path
#' @param name name
#' @author Angela Andreella
#' @return lambda
#' @export
#' @importFrom stats pbeta
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @importFrom graphics legend

lambdaPlot <- function(copes = NULL, family, ct = c(0,1), alpha, delta = NULL, P = NULL, mask, path = getwd(), name = "plot"){
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  if(is.null(copes) & is.null(P)){stop('Please insert pvalues matrix or copes images')}
  if(!is.null(P) & is.unsorted(P[1,])){pvalues <- rowSortC(P)}
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
    res <- signTest(X=scores, B = 1000,alternative = "two.sided", rand = F) #variables times number of permutation
    
    pvalues <- cbind(res$pv,res$pv_H0)
    pvalues = t(pvalues)
    rm(res)
    rm(scores)
    rm(copes)
    rm(img)
    
    pvalues <- rowSortC(pvalues)
    
    
  }
  
  
  l <- c()
  w <- dim(pvalues)[1]
  m <- dim(pvalues)[2]
  
  if(is.null(delta) ){delta = 0}
  for(j in 1:w){
    minc <- sum(pvalues[j,] <=min(ct)) + 1
    maxc <- sum(pvalues[j,] <=max(ct))
    if(family =="simes"){
      
      minc = minc + delta
      lambda <- ((m-delta)*(pvalues[j,minc:maxc]))/((c(minc:maxc)-delta)*alpha)
    }
    if(family == "beta"){
      lambda <- pbeta(q =pvalues[j,c(minc:maxc)],shape1 = c(minc:maxc),shape2 =m+1-c(minc:maxc))
      
    }
    
    if(family == "finner"){
      minc = minc + delta
      lambda <- (pvalues[j,c(minc:maxc)]*(m - 1) )/ (alpha * (c(minc:maxc) - delta) * (1 - pvalues[j,c(minc:maxc)]))
      
    }
    
    
    if(family =="higher.criticism"){
      
      lambda <- (sqrt(m)*((c(minc:maxc)/m) - pvalues[j,c(minc:maxc)]))/(sqrt(pvalues[j,c(minc:maxc)]*(1-pvalues[j,c(minc:maxc)])))
      
    }
    if(family=="beta"){
      l[j] <- min(lambda[lambda>0])
    }else{
      l[j] <- min(lambda)
    }
    
  }
  
  lambdaE <- sort(l)[floor(alpha*w)+1]
  cvE<- cv(pvalues = pvalues, family = family, alpha = alpha, lambda = lambdaE, delta = delta)
  lcvL <- function(family,delta=NULL, lambda){
    cvO<- cv(pvalues = pvalues, family = family, alpha = alpha, lambda = lambda, delta = delta)
    lines(cvO, lwd =1, col= "grey")
  }
  lcvLV <- Vectorize(lcvL,vectorize.args = c("lambda"))
  
  png(paste0(path,"/", name, ".png")) 
  plot(pvalues[1,], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]))
  lcvLV(family = family, delta = delta, lambda = l)
  for(i in 2:nrow(pvalues)){
    
    lines(pvalues[i,],col='black',type="l")
    
  }
  lines(pvalues[1,], lwd =2, col= 'red')
  lines(cvE, col= 'red', lwd =2, lty = "dashed")
  legend('top',c("Observed Pvalues", "Null Pvalues", "Optimal Critical Vector", "Critical Vectors"), col= c("red", "black", "red", "grey"),lty=c(1,1,2,1),lwd =2)

  dev.off()
  


}


#lambdaOptAprox <- function(pvalues, family, ct = tS, alpha, shift = delta,cb, Kc){
#  l <- c()
#  w <- dim(pvalues)[1]
#  m <- dim(pvalues)[2]
#  quant <- (pvalues+shift)*m/(alpha)
#  lambdaA <- c()
#  for(j in 1:w){
#    minc <- sum(pvalues[j,Kc[,cb]] <=min(ct)) + 1
#    maxc <- sum(pvalues[j,Kc[,cb]] <=max(ct))
#    
#    if(is.null(shift)){shift = 0}
#    lambdaA[j] <- min(sapply(c(minc:maxc), function(x) quant[j, Kc[   sort.list(pvalues[j,Kc[,cb]])[x] ,cb   ] ]/x  ))
#    #lamb <- min(lamb,  betaquantsS[, combs2[   sort.list(pvmatr.uns[j,combs2[,c]]),c]][j,a]/a  )


#  }


# nRej_Perm <- rejPerm(pvalues,min(ct))
# quantRej_min <- sort(nRej_Perm)[ceiling((1-alpha)*w)]
# if(quantRej_min>0){
#   lambdaE <- min(lambdaE, (min(ct) + shift)*m/(quantRej_min*alpha))}

#  return(lambdaA)
#}

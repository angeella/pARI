#' @title ARI Permutation-based 
#' @description Performs ARI using permutation local test
#' @param copes = list of copes map
#' @param thr = threshold to construct cluster map
#' @param mask = mask map, niftiImage class object or path
#' @param alpha = alpha level
#' @param delta = sdo you want to consider at least delta size set?
#' @param summary_stat = Choose among \code{=c("max", "center-of-mass")}
#' @param silent \code{FALSE} by default.
#' @family = if permutation approach, which family for the confidence envelope?
#' @B = number of permutation, default 1000
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{discoveries} number of discoveries in the set selected, cluster id, maximum test statistic and relative coordinates
#' @export

ARIpermCT <- function(copes, thr, mask=NULL, alpha=.1, clusters = NULL,
                      summary_stat=c("max", "center-of-mass"),silent=FALSE, family = NULL, delta = NULL, B = 1000, ...){
  
  if(is.character(mask)){mask = RNifti::readNifti(mask)}
  if(!is.list(copes)){stop("Please insert the list of copes as list class object")}

  img_dims <- c(91,  109 , 91)
  img <- array(NA, c(img_dims, length(copes)))
  
  for (sid in 1:length(copes)) {  
    img[,,,sid] <- copes[[sid]]
    
  }
  
  scores <- matrix(img,nrow=(91*109*91),ncol=length(copes))
  resO <-oneSample(X=scores,alternative = "two.sided")
  
  scores <- scores[which(mask==1),]
  res <- signTest(X=scores, B = B,alternative = "two.sided") #variables times number of permutation
  
  pvalues <- cbind(res$pv,res$pv_H0)
  pvalues = t(pvalues)
  Statmap = array(data = resO$Test, dim = c(91,109,91))
  Statmap[!mask]=0
  rm(res)
  rm(scores)
  rm(copes)
  rm(img)
  
  if(is.null(clusters) & !is.null(thr)){clusters <- cluster_threshold(Statmap>thr)} 
  #clusters = get_array(clusters,map_dims=dim(Pmap))

  # called=match.call()
  summary_stat=match.arg(summary_stat,c("max", "center-of-mass"))
  
  # get the indices of the mask
  mask=which(mask!=0) #vector of n voxels
  
  #As first, we compute the optimal lambda
  #pv_id <- Pmap[mask]
  #pvalues <- pvalues[,mask]
  #pvalues <- rbind(pv_id,pvalues)
  #pvalues_ord <- t(apply(pvalues, 1, sort))
    
  pvalues_ord <- rowSortC(pvalues)
  praw <- pvalues_ord[1,]
  lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, delta = delta) 
  #cvh <- cvhPerm(praw = praw, alpha = alpha, shift = shift, family = family, lambda = lambda)
  #cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
  #cv <- sapply(c(1:length(praw)), function(x) (((x-8) * alpha * lambda)/length(praw)))
  cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda)
 
  # define number of clusters
  clstr_id=sort(unique(as.vector(clusters[mask])),decreasing = TRUE)
  
  if(is.function(Statmap)) {
    StatFun=Statmap
  } else {
    #Statmap= get_array(Statmap,map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
  }
  
  #apply summaries to each cluster (and all the rest in an extra cluster)
    out=plyr::laply(clstr_id,function(i){
      ix=clusters==i
      ix[-mask]=FALSE
      
      cluster_ids=which(ix,arr.ind = TRUE)
      cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
      #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
      #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
      #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
      unlist(c(summary_perm_roi(cv = cvOpt,ix=ix[mask],pvalues = pvalues),
               summary_cluster(cluster_ids)[-1])
      )
    })
  rownames(out)=paste("cl",sep="",clstr_id)
  
  # attr(out,"call")=called
  if(!silent) print(out)
  out
}














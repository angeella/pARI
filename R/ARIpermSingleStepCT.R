#' @title ARI Permutation-based 
#' @description Performs ARI using permutation local test
#' @usage ARIpermCT(copes, thr, mask, alpha, clusters,summary_stat,silent, family, delta, B, ct)
#' @param copes list of copes map
#' @param thr threshold to construct cluster map
#' @param mask mask map, niftiImage class object or path
#' @param alpha alpha level
#' @param clusters clusters map as niftiImage class object or path, if NULL it is computed considering the threshold 3.2
#' @param summary_stat Choose among "max", "center-of-mass"
#' @param silent FALSE by default.
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is simes
#' @param delta do you want to consider at least delta size set?
#' @param B number of permutation, default 1000
#' @param ct set of thresholds
#' @param rand logical. Should p values computed by permutation distribution?
#' @author Angela Andreella
#' @return Returns a list with the following objects: discoveries number of discoveries in the set selected, cluster id, maximum test statistic and relative coordinates
#' @export
#' @importFrom RNifti readNifti
#' @importFrom plyr laply

ARIpermCT <- function(copes, thr=NULL, mask=NULL, alpha=.1, clusters = NULL,
                      summary_stat=c("max", "center-of-mass"),silent=FALSE, family = "simes", delta = NULL, B = 1000, ct = c(0,1), rand = FALSE){
  
  "%ni%" <- Negate("%in%")
  #check alpha
  val_alpha = sapply(c(1:B), function(x) (B-x)/B)
  if(alpha %ni% val_alpha){stop('please insert valid values for alpha and B')}
  
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  
  family <- match.arg(tolower(family), family_set)
  if(is.character(mask)){mask = readNifti(mask)}
  if(!is.list(copes)){stop("Please insert the list of copes as list class object")}

  img_dims <- c(91,  109 , 91)
  img <- array(NA, c(img_dims, length(copes)))
  
  for (sid in 1:length(copes)) {  
    img[,,,sid] <- copes[[sid]]
    
  }
  
  scores <- matrix(img,nrow=(91*109*91),ncol=length(copes))
  scores[!mask,] = NA
  resO <-oneSample(X=scores,alternative = "two.sided")
  
  scores <- scores[which(mask==1),]
  res <- signTest(X=scores, B = B,alternative = "two.sided", rand = rand) #variables times number of permutation
  
  pvalues <- cbind(res$pv,res$pv_H0)
  pvalues = t(pvalues)
  Statmap = array(data = resO$Test, dim = c(91,109,91))
  Statmap[!mask]=0
  rm(res)
  rm(scores)
  rm(copes)
  rm(img)
  
  if(is.null(clusters) & !is.null(thr)){clusters <- cluster_threshold(Statmap>thr)}
  if(!is.null(clusters) & is.null(thr)){
    if(is.character(clusters)){
      clusters = readNifti(clusters)
    }else{
      clusters = get_array(clusters)
    }

      clusters = array(clusters,dim(clusters))
  } 
  if(is.null(clusters) & is.null(thr) & !is.null(mask)){clusters <- array(mask,dim(mask))}
  if(is.null(clusters) & is.null(thr) & is.null(mask)){stop("Please insert mask, threshold value or cluster map")}
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
  cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda, delta = delta)
 
  # define number of clusters
  clstr_id=sort(unique(as.vector(clusters[mask])),decreasing = TRUE)
  
  if(is.function(Statmap)) {
    StatFun=Statmap
  } else {
    #Statmap= get_array(Statmap,map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
  }
  
  #apply summaries to each cluster (and all the rest in an extra cluster)
    out=laply(clstr_id,function(i){
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
  if(!is.null(dim(out))){
    rownames(out)=paste("cl",sep="",clstr_id)
  }
  
  
  # attr(out,"call")=called
  if(!silent) print(out)
  return(list(out = out,clusters = clusters))
}














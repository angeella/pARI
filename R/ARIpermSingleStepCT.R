#' @title ARI Permutation-based 
#' @description Performs ARI using permutation local test
#' @param Pmap = pvalues map
#' @param thr = threshold to construct cluster map
#' @param mask = mask mmap
#' @param alpha = alpha level
#' @param delta = sdo you want to consider at least delta size set?
#' @param Statmap = test statistic map
#' @summary_stat = Choose among \code{=c("max", "center-of-mass")}
#' @param silent \code{FALSE} by default.
#' @type = parametric o permutation approach?
#' @family = if permutation approach, which family for the confidence envelope?
#' @pvalues = if permutation approach, which family for the confidence envelope?
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{discoveries} number of discoveries in the set selected, \code{pvalues} raw pvalues
#' @export

ARIpermCT <- function(Pmap, thr, mask=NULL, alpha=.1,Statmap=function(ix) -qnorm(Pmap[ix]),
                      summary_stat=c("max", "center-of-mass"),silent=FALSE, type="perm", family = NULL, pvalues = NULL, delta = NULL, ...){
  
  # get, fix, check parameters inconsistencies
  Pmap = get_array(Pmap)
  Statmap = get_array(Statmap)
  mask = get_array(mask,map_dims=dim(Pmap))
  Statmap[!mask]=0
  clusters <- cluster_threshold(Statmap>thr)
  #clusters = get_array(clusters,map_dims=dim(Pmap))
  if(!is.null(pvalues)){
    if(is.character(pvalues)){load(paste0(pvalues))}else{
      if(!is.matrix(pvalues))
        stop("Please insert path or matrix pvalues perm")
        }
    
    pvalues = t(pvalues)}
  if(is.function(Statmap)) {
    StatFun=Statmap
  } else {
    Statmap= get_array(Statmap,map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
  }
  
  # called=match.call()
  summary_stat=match.arg(summary_stat,c("max", "center-of-mass"))
  
  # get the indices of the mask
  mask=which(mask!=0) #vector of n voxels
  
  if(type=='perm'){
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
  }else{
    #perform hommel
    hom <-hommel::hommel(Pmap[mask])
    if(!silent) {temp=(summary(hom))
    cat("\n")}
  }
  
  
  # define number of clusters
  clstr_id=sort(unique(as.vector(clusters[mask])),decreasing = TRUE)
  
  #apply summaries to each cluster (and all the rest in an extra cluster)
  if(type=='perm'){
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
  }else{
    out=plyr::laply(clstr_id,function(i){
      ix=clusters==i
      ix[-mask]=FALSE
      
      cluster_ids=which(ix,arr.ind = TRUE)
      cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
      
      unlist(c(summary_hommel_roi(hommel = hom,ix=ix[mask], alpha= alpha),
               summary_cluster(cluster_ids)[-1])
      )
    })}
  rownames(out)=paste("cl",sep="",clstr_id)
  
  # attr(out,"call")=called
  if(!silent) print(out)
  out
}














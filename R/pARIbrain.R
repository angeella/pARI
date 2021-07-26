#' @title Permutation-based All-Resolutions Inference for brain imaging
#' @description The main function for single step All-Resolutions Inference (ARI) method based on critical vectors constructed by permutations for fMRI cluster analysis. 
#' @usage pARIbrain(copes, thr, mask, alpha, clusters, 
#' alternative, summary_stat, silent, family, delta, B, rand,
#' step.down, max.step)
#' @param copes The list of copes, i.e., constrasts maps, one for each subject used to compute the statistical tests.
#' @param thr threshold used to construct the cluster map.
#' @param mask 3D array of locicals (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @param alpha alpha level
#' @param clusters 3D array of cluster ids (0 when voxel does not belong to any cluster) or a (character) nifti file name. 
#' If \code{cluster=NULL} the cluster map is computed by the \code{\link{cluster_threshold}} function with threshold equals 3.2.
#' @param alternative a character string referring to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param summary_stat Choose among \code{=c("max", "center-of-mass")}. 
#' @param silent \code{FALSE} by default.
#' @param family by default \code{family="simes"}. Choose a family of confidence envelopes to compute the critical vector from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.
#' @param delta by default \code{delta = 0}. Do you want to consider sets with at least delta size?
#' @param B by default \code{B = 1000}. Number of permutations.
#' @param rand by default \code{rand = FALSE}. 
#' A logical value, if \code{rand = TRUE}, the pvalues are computed by \code{\link{rowRanks}}.
#' @param iterative if \code{iterative = TRUE}, the iterative iterative method for improvement of confidence envelopes is applied. Default is \code{FALSE}.
#' @param approx if \code{iterative = TRUE} and you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb if \code{approx = TRUE}, you must decide how many large random subcollection (level of approximation) considered.
#' @param step.down by default \code{step.down = FALSE}. If you want to compute the lambda calibration parameter using the step down approach put TRUE.
#' @param max.step by default \code{max.step = 10}. Maximum number of steps for the step down approach
#' @author Angela Andreella
#' @return Returns a matrix with the following objects: 
#' Size, FalseNull, TrueNull, ActiveProp and other statistics for each cluster.
#' @export
#' @importFrom RNifti readNifti
#' @importFrom plyr laply
#' @importFrom ARIbrain cluster_threshold
#' @references For the general framework of All-Resolutions Inference see:
#' 
#' Goeman, Jelle J., and Aldo Solari. "Multiple testing for exploratory research." Statistical Science 26.4 (2011): 584-597.
#'
#' For All-Resolutions Inference for functional Magnetic Resonance Imaging data see: 
#' 
#' Rosenblatt, Jonathan D., et al. "All-resolutions inference for brain imaging." Neuroimage 181 (2018): 786-796.
#' 
#'
#' For permutation-based All-Resolutions Inference see:
#' 
#' Andreella, Angela, et al. "Permutation-based true discovery proportions for fMRI cluster analysis." arXiv preprint arXiv:2012.00368 (2020).
#' 
#' @examples
#' \dontrun{
#' library(fMRIdata)
#' data(Auditory_clusterTH3_2)
#' data(Auditory_copes)
#' data(Auditory_mask)
#' auditory_out <- pARIbrain(copes = Auditory_copes, 
#' cluster = Auditory_clusterTH3_2, mask = Auditory_mask, 
#' alpha = 0.05, silent = TRUE)
#' auditory_out$out
#' }


pARIbrain <- function(copes, thr=NULL, mask=NULL, alpha=.05, clusters = NULL, 
                      alternative = "two.sided", summary_stat=c("max", "center-of-mass"),
                      silent=F, family = "simes", delta = 0, B = 1000, rand = FALSE, 
                      iterative = FALSE, approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10){
  
  "%ni%" <- Negate("%in%")
  #check alpha
  val_alpha = sapply(c(1:B), function(x) (B-x)/B)
  if(!(alpha %in% val_alpha)){stop('please insert valid values for alpha and B')}
  
  family_set <- c("simes", "aorc", "beta", "higher.criticism")
  alternative_set <- c("two.sided", "greater", "lower")
  family <- match.arg(tolower(family), family_set)
  alternative <- match.arg(tolower(alternative), alternative_set)
  
  if(is.character(mask)){mask = readNifti(mask)}
  if(!is.list(copes)){stop("Please insert the list of copes as list class object")}

  img_dims <- c(91,  109 , 91)
  img <- array(NA, c(img_dims, length(copes)))
  
  for (sid in 1:length(copes)) {  
    img[,,,sid] <- copes[[sid]]
    
  }
  
  scores <- matrix(img,nrow=(91*109*91),ncol=length(copes))
  scores[!mask,] = NA
  resO <-oneSample(X=scores,alternative = alternative)
  
  scores <- scores[which(mask==1),]
  res <- signTest(X=scores, B = B,alternative = alternative, rand = rand) #variables times number of permutation
  
  pvalues <- cbind(res$pv,res$pv_H0)
#  pvalues = t(pvalues)
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
    
#  pvalues_ord <- rowSortC(pvalues)
#  praw <- pvalues_ord[1,]

  lambda <- lambdaOpt(pvalues = pvalues, family = family, alpha = alpha, delta = delta, step.down = step.down, max.step = max.step) 
  #cvh <- cvhPerm(praw = praw, alpha = alpha, shift = shift, family = family, lambda = lambda)
  #cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
  #cv <- sapply(c(1:length(praw)), function(x) (((x-8) * alpha * lambda)/length(praw)))
  cvOpt = criticalVector(pvalues = pvalues, family = family, alpha = alpha, lambda= lambda, delta = delta)
 
  # define number of clusters
  clstr_id=sort(unique(as.vector(clusters[mask])),decreasing = TRUE)
  
  if(is.function(Statmap)) {
    StatFun=Statmap
  } else {
    #Statmap= get_array(Statmap,map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
  }
  clstr_id <- clstr_id[clstr_id!=0]
  #apply summaries to each cluster (and all the rest in an extra cluster)
    out=laply(clstr_id,function(i){
      print(i)
      ix=clusters==i
      ix[-mask]=FALSE
      
      cluster_ids=which(ix,arr.ind = TRUE)
      cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
      #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
      #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
      #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
      #print(dim(pvalues))
      unlist(c(summary_perm_roi(cv = cvOpt,ix=ix[mask],pvalues = pvalues, 
                                iterative = iterative, approx = approx, ncomb = ncomb, 
                                family = family, delta = delta, alpha = alpha),
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














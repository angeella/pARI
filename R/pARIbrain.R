#' @title Permutation-based All-Resolutions Inference for brain imaging.
#' @description The main function for brain imaging All-Resolutions Inference (ARI) method based on critical vectors constructed 
#' using the p-values permutation distribution. The function computes simultaneous lower bounds for the number of true discoveries 
#' for each set of hypotheses specified in \code{ix} controlling family-wise error rate.
#' @usage pARIbrain(copes, thr=NULL, mask=NULL, alpha=.05, clusters = NULL, 
#' alternative = "two.sided", summary_stat=c("max", "center-of-mass"),
#' silent=FALSE, family = "simes", delta = 0, B = 1000, rand = FALSE, 
#' iterative = FALSE, approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...)
#' @param copes list of NIfTI file. The list of copes, i.e., constrasts maps, one for each subject used to compute the statistical tests.
#' @param thr numeric value. Threshold used to construct the cluster map. Default @NULL.
#' @param mask NIfTI file or character string. 3D array of logical values (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @param alpha numeric value in `[0,1]`. It expresses the alpha level to control the family-wise error rate. Default 0.05.
#' @param clusters NIfTI file or character string. 3D array of cluster ids (0 when voxel does not belong to any cluster) or a (character) NIfTI file name. 
#' If \code{cluster=NULL} the cluster map is computed by the \code{\link{cluster_threshold}} function with threshold equals \code{thr}.
#' @param alternative character string. It refers to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param summary_stat character string. Choose among \code{=c("max", "center-of-mass")}. 
#' @param silent Boolean value. Default @FALSE. If @TRUE the function prints the results.
#' @param family string character. Choose a family of confidence envelopes to compute the critical vector 
#' from \code{"simes"}, \code{"aorc"}, \code{"beta"} and \code{"higher.criticism"}.#' @param alpha alpha level.
#' @param delta numeric value. It expresses the delta value, please see the references. Default to 0. 
#' @param B numeric value. Number of permutations, default to 1000. 
#' @param rand Boolean value. Default @FALSE. If \code{rand = TRUE}, the p-values are computed by \code{\link{rowRanks}}.
#' @param iterative Boolean value. If \code{iterative = TRUE}, the iterative method for improvement of confidence envelopes is applied. Default @FALSE.
#' @param approx Boolean value. Default @TRUE. If you are treating high dimensional data, we suggest to put \code{approx = TRUE} to speed up the computation time.
#' @param ncomb Numeric value. If \code{approx = TRUE}, you must decide how many random subcollections (level of approximation) considered.
#' @param step.down Boolean value. Default @FALSE If you want to compute the lambda calibration parameter using the step-down approach put @TRUE.
#' @param max.step Numeric value. Default to 10. Maximum number of steps for the step down approach, so useful when \code{step.down = TRUE}.
#' @param ... further arguments. See \code{signTest}.
#' @author Angela Andreella
#' @return A list with elements
#' - out: data.frame containing the size, the number of false null hypotheses, 
#' the number of true null hypotheses, the lower bound for the true discovery proportion, and other
#' statistics for each cluster.
#' - clusters: matrix describing the clusters analyzed.
#' 
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
#' For permutation-based All-Resolutions Inference see:
#' 
#' Andreella, Angela, et al. "Permutation-based true discovery proportions for fMRI cluster analysis." arXiv preprint arXiv:2012.00368 (2020).
#' 
#' @examples
#' \dontrun{
#' library(remotes)
#' install_github("angeella/fMRIdata")
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
                      silent=FALSE, family = "simes", delta = 0, B = 1000, rand = FALSE, 
                      iterative = FALSE, approx = TRUE, ncomb = 100, step.down = FALSE, max.step = 10, ...){
  
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
  resO <-oneSamplePar(X=scores,alternative = alternative)
  
  scores <- scores[which(mask==1),]
  res <- signTest(X=scores, B = B,alternative = alternative, rand = rand, ...) #variables times number of permutation
  
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

  lambda <- lambdaOpt(pvalues = pvalues, family = family, alpha = alpha, delta = delta, step.down = step.down, max.step = max.step) 
  if(lambda == 0){lambda <- 0.05}
  #and critical vector
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














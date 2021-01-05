#' @title Create Statistical Parametric Mapping (SPM)
#' @description \code{Statmap} is used to create the statistical parametric mapping in NIfTI format.
#' @usage Statmap(copes,alternative,path,name, Pmap, mask)
#' @param copes The list of copes, i.e., constrasts maps, one for each subject used to compute the statistical tests.
#' @param alternative a character string referring to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param path path used to save the NIfTI file, the path does not must end with \code{/}.
#' @param name choose the name of your NIfTI file
#' @param Pmap by default \code{Pmap=FALSE}. If \code{TRUE} the SPM of the pvalues is returned.
#' @param mask 3D array of locicals (i.e. \code{TRUE/FALSE} in/out of the brain). Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @author Angela Andreella
#' @return the statistical parametric mapping in NIfTI format is saved in the path specified, if \code{path = NULL} the current working directory of the \code{R} process is used.
#' @export
#' @author Angela Andreella
#' @importFrom RNifti writeNifti
#' @importFrom RNifti readNifti
#' @examples
#' \dontrun{
#' library(fMRIdata)
#' data(Auditory_copes)
#' data(Auditory_mask)
#' Statmap(copes = Auditory_copes, mask = Auditory_mask)
#' }

Statmap <- function(copes, alternative = "two.sided", path = getwd(), name = "map", Pmap = FALSE, mask = NULL){
  
  if(!is.list(copes)){stop("Please insert the list of copes as list class object")}
  if(is.character(mask)){mask = readNifti(mask)}
  
  img_dims <- c(91,  109 , 91)
  img <- array(NA, c(img_dims, length(copes)))
  
  for (sid in 1:length(copes)) {  
    img[,,,sid] <- copes[[sid]]
    
  }
  
  scores <- matrix(img,nrow=(91*109*91),ncol=length(copes))
  scores[!mask,] = NA
  resO <-oneSample(X=scores,alternative = alternative)
  if(Pmap){
    pv = array(data = resO$p, dim = c(91,109,91))
    writeNifti(pv,file = paste0(path, "/P", name,".nii.gz"),template = mask)
  }

  t = array(data = resO$T, dim = c(91,109,91))
  writeNifti(t,file = paste0(path, "/Stat", name,".nii.gz"),template = mask)
}


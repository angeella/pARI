#' @title Create Statistical Parametric Mapping (SPM)
#' @description It creates the statistical parametric mapping in NIfTI format.
#' @usage Statmap(copes, alternative = "two.sided", path = getwd(), 
#' name = "map", Pmap = FALSE, mask = NULL)
#' @param copes list of NIfTI file. The list of copes, i.e., constrasts maps, one for each subject used to compute the statistical tests.
#' @param alternative character string. It refers to the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"lower"}.
#' @param path character string. Path to save the NIfTI file. The path does not must end with \code{/}.
#' @param name character string. The name of the map NIfTI file that will be used.
#' @param Pmap Boolean value. If \code{TRUE} the SPM of the pvalues is returned. Default @FALSE.
#' @param mask NIfTI file or character string. 3D array of logical values (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @author Angela Andreella
#' @return Save the Statistical Parametric Mapping Nifti file in \code{path} with name specified in \code{name}.
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
  
  #check copes
  if(!is.list(copes)){stop("Please insert the list of copes as list class object")}

  img_dims <- dim(copes[[1]])
  
  #check mask
  if(!is.null(mask)){
    if(!is.character(mask) && !is.array(mask)){stop("mask must be an array or a path")}
    if(is.character(mask)){mask = readNifti(mask)}
    if(!all(dim(mask) == img_dims)){stop("incompatible dimensions of mask and copes")}
  }else{
    mask <- array(1, img_dims)
  }
  
  #create score matrix
  img <- array(NA, c(img_dims, length(copes)))
  for (sid in 1:length(copes)) {  
    img[,,,sid] <- copes[[sid]]
    
  }
  scores <- matrix(img,nrow=(img_dims[1]*img_dims[2]*img_dims[3]),ncol=length(copes))
  scores[!mask,] = NA
  
  resO <-oneSamplePar(X=scores,alternative = alternative)
  if(Pmap){
    pv = array(data = resO$p, dim = img_dims)
    writeNifti(pv,file = paste0(path, "/P", name,".nii.gz"),template = mask)
  }

  t = array(data = resO$T, dim = img_dims)
  writeNifti(t,file = paste0(path, "/Stat", name,".nii.gz"),template = mask)
}


#' @title Statmap
#' @description create Statmap nifti image
#' @usage Statmap(copes,alternative,path,name, Pmap, mask)
#' @param copes list of copes
#' @param alternative which alternative hypothesis perform, default is two.sided
#' @param path path to save the nifti file, the path doesn' t must end with /
#' @param name choose the name of your nifti file
#' @param Pmap if you want also the Pvalue map insert TRUE, default is FALSE
#' @param mask mask map
#' @author Angela Andreella
#' @return lambda
#' @export
#' @importFrom RNifti writeNifti
#' @importFrom RNifti readNifti


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


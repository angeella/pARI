#' @title True Discovery Proportion brain map
#' @description Performs True Discovery Proportion brain map. 
#' @usage map_TDP(ARIout,path,name,mask)
#' @param ARIout output object by \code{pARIbrain}.
#' @param path path used to save the NIfTI file, the path does not must end with \code{/}.
#' @param name choose the name of the NIfTI file.
#' @param mask 3D array of locicals (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @author Angela Andreella
#' @return Returns a TDP map
#' @export
#' @importFrom RNifti writeNifti
#' @importFrom RNifti readNifti

map_TDP <- function(ARIout,path= getwd(), name = "tdp", mask){
  
  
  if(is.character(mask)){mask = readNifti(mask)}
  
  cl = as.numeric(gsub("cl", "",names(ARIout$out[,4])))
  tdp = as.numeric(ARIout$out[,4])
  tdp[length(cl)] = 0
  
  map = ARIout$clusters
  
  if(!is.array(map)){ map <- array(map,dim(map))}
    
  for(i in 1:length(cl)){
    
    map[map==cl[i]]=tdp[i]
  }
  
  writeNifti(map,file = paste0(path, "/", name, ".nii"),template = mask)
  
}






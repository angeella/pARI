#' @title True Discovery Proportion brain map
#' @description Performs the True Discovery Proportion brain map. 
#' @usage map_TDP(ARIout,path,name,mask)
#' @param ARIout output object by \code{\link{pARIbrain}}.
#' @param path character string, path to save the NIfTI file. The path does not must end with \code{/}.
#' @param name character string, the name of the map NIfTI file that will be used.
#' @param mask 3D array of logicals (i.e. \code{TRUE/FALSE} in/out of the brain). 
#' Alternatively it may be a (character) NIfTI file name. If \code{mask=NULL}, 
#' it is assumed that none of the voxels have to be excluded.
#' @author Angela Andreella
#' @return Returns the True Discovery Proportion NIfTI map.
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






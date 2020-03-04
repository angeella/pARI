#' @title map True Discovery Proportion
#' @description Performs brain map True Discovery Proportion
#' @usage map_TDP(ARIout,path,name,mask)
#' @param ARIout output object ARI permutation
#' @param path path
#' @param name name
#' @param mask mask map, niftiImage class object or path
#' @author Angela Andreella
#' @return Returns a TDP map
#' @export
#' @importFrom RNifti writeNifti

map_TDP <- function(ARIout,path= getwd(), name = "tdp", mask){
  
  if(is.character(mask)){mask = RNifti::readNifti(mask)}
  
  cl = as.numeric(gsub("cl", "",names(ARIout$out[,4])))
  tdp = as.numeric(ARIout$out[,4])
  tdp[length(cl)] = 0
  
  map = ARIout$clusters
  
  for(i in 1:length(cl)){
    
    map[map==cl[i]]=tdp[i]
  }
  
  writeNifti(map,file = paste0(path, "/", name, ".nii"),template = mask)
  
}






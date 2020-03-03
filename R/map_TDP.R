#' @title map True Discovery Proportion
#' @description Performs brain map True Discovery Proportion
#' @param ARIout = output object ARI permutation
#' @param B = number of permutations to perform, default is 1000.
#' @param alternative = character referring to the alternative hypothesis, "two.sided", "greater" or "less". Default is "two.sided"
#' @param seed = specify seed, default is 1234.
#' @param mask = mask map, niftiImage class object or path
#' @author Angela Andreella
#' @return Returns a list with the following objects: \code{Test} observed one sample t-test, \code{Test_H0} Test statistics under H0, \code{pv} observed p-values, \code{pv_H0} p-values under H0
#' @export

map_TDP <- function(ARIout,path= getwd(), name = "tdp", mask){
  
  if(is.character(mask)){mask = RNifti::readNifti(mask)}
  
  cl = as.numeric(gsub("cl", "",names(ARIout$out[,4])))
  tdp = as.numeric(ARIout$out[,4])
  tdp[length(cl)] = 0
  
  map = ARIout$clusters
  
  for(i in 1:length(cl)){
    
    map[map==cl[i]]=tdp[i]
  }
  
  RNifti::writeNifti(map,file = paste0(path, "/", name, ".nii"),template = mask)
  
}






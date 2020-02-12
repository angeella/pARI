summary_cluster <- function(coord_and_values,summary_stat=c("max", "center-of-mass")){
  # compute max and/or centre of gravity, see below
  #no_Inf = which(as.vector(coord_and_values[,4]) != Inf)
  #coord_and_values = as.data.frame(t(coord_and_values[no_Inf,]))
  name_stat=names(coord_and_values)[1]
  summary_stat=match.arg(summary_stat,c("max", "center-of-mass"))
  out=list(Size=nrow(coord_and_values))
  if(summary_stat=="max"){
    id_max=which.max(coord_and_values[,4])
    out=c(out,coord_and_values[id_max,])
  } else if(summary_stat=="center-of-mass"){
    id_mean=colMeans(coord_and_values[,-1,drop=FALSE])
    id_closest_to_baricenter=which.min(rowSums(t(t(coord_and_values[,-1,drop=FALSE])-id_mean)^2))
    out=c(out,coord_and_values[id_closest_to_baricenter,])
    names(out)[names(out)==name_stat]=summary_stat
  }
  out
}
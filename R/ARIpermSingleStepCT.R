###################################ARIbrain SingleStepCT#####################################################################


ARIpermCT <- function(Pmap, clusters, mask=NULL, alpha=.1,Statmap=function(ix) -qnorm(Pmap[ix]),
                      summary_stat=c("max", "center-of-mass"),silent=FALSE, type="perm", family = NULL, pvalues = NULL, shift = NULL, ...){

  # get, fix, check parameters inconsistencies
  Pmap = get_array(Pmap)
  clusters = get_array(clusters,map_dims=dim(Pmap))
  mask = get_array(mask,map_dims=dim(Pmap))
  if(is.function(Statmap)) {
    StatFun=Statmap
  } else {
    Statmap= get_array(Statmap,map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
  }
  
  # called=match.call()
  summary_stat=match.arg(summary_stat,c("max", "center-of-mass"))
  
  # get the indices of the mask
  mask=which(mask!=0) #vector of n voxels
  
  if(type=='perm'){
    #As first, we compute the optimal lambda
    #pv_id <- Pmap[mask]
    #pvalues <- pvalues[,mask]
    #pvalues <- rbind(pv_id,pvalues)
    pvalues_ord <- t(apply(pvalues, 1, sort))
    
    #pvalues_ord <- rowSortC(pvalues)
    praw <- pvalues_ord[1,]
    lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, shift = shift) 
    #cvh <- cvhPerm(praw = praw, alpha = alpha, shift = shift, family = family, lambda = lambda)
    #cv <- sapply(c(1:length(praw)), function(x) ((x * alpha * lambda)/length(praw))- shift)
    #cv <- sapply(c(1:length(praw)), function(x) (((x-8) * alpha * lambda)/length(praw)))
    cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, shift = shift, lambda= lambda)
  }else{
    #perform hommel
    hom <-hommel::hommel(Pmap[mask])
    if(!silent) {temp=(summary(hom))
    cat("\n")}
  }
  
  
  # define number of clusters
  clstr_id=sort(unique(as.vector(clusters[mask])),decreasing = TRUE)
  
  #apply summaries to each cluster (and all the rest in an extra cluster)
  if(type=='perm'){
    out=plyr::laply(clstr_id,function(i){
      ix=clusters==i
      ix[-mask]=FALSE
      
      cluster_ids=which(ix,arr.ind = TRUE)
      cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
      #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
      #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
      #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
      unlist(c(summary_perm_roi(cv = cvOpt,ix=ix[mask]),
               summary_cluster(cluster_ids)[-1])
      )
    })
  }else{
    out=plyr::laply(clstr_id,function(i){
      ix=clusters==i
      ix[-mask]=FALSE
      
      cluster_ids=which(ix,arr.ind = TRUE)
      cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
      
      unlist(c(summary_hommel_roi(hommel = hom,ix=ix[mask], alpha= alpha),
               summary_cluster(cluster_ids)[-1])
      )
    })}
  rownames(out)=paste("cl",sep="",clstr_id)
  
  # attr(out,"call")=called
  if(!silent) print(out)
  out
}

summary_hommel_roi <- function(hommel,ix,alpha){
  Total=length(hommel@p[ix])
  False_Null=hommel::discoveries(hommel, alpha=alpha, ix=ix)
  True_Null=Total-False_Null
  Active_Proportion= tdp(hommel, ix=ix, alpha = alpha)
  list(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
}

summary_perm_roi <- function(cv,ix){
  idix <- which(ix)
  p <- pvalues[1,idix]
  Total = length(p)
  False_Null= dI(ix = idix,cv = cv,praw = pvalues[1,])
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  list(Size=Total,FalseNull=False_Null,TrueNull=True_Null,ActiveProp=Active_Proportion)
}



tdp <- function (hommel, ix, alpha = alpha) 
{
  m <- length(hommel@p)
  if (missing(ix)) {
    d <- discoveries(hommel, alpha = alpha)
    k <- m
  }
  else {
    p <- hommel@p[ix]
    k <- length(p)
    d <- discoveries(hommel, ix, alpha = alpha)
  }
  d/k
}


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


get_array <- function(map,map_dims=NULL){
  if(is.null(map)){
    if(is.null(map_dims)) {
      stop("The dims of map are not defined")
    } else {
      map=array(TRUE,map_dims)
      return(map)
    }
  }
  
  if(is.character(map))
    map=RNifti::readNifti(map)
  
  if(!is.null(map_dims)) if(any(dim(map)!=map_dims)) stop("The dims of map: ",
                                                          paste(dim(map),sep="x"),
                                                          " don't fit the dims of map_dims: ",paste(map_dims,sep="x"))
  map
}

cluster_threshold <- function(map, max_dist=sqrt(3)){
  ### slower:
  # map=spmT
  # threshold=3.2
  # nmat <- expand.grid(-1:1, -1:1, -1:1)
  # nmat <- nmat[-c(1,3,7,9,14,19,21,25,27), ]
  # system.time(
  # {Suprathreshold_TF = cluster.threshold(spmT>=3.2, nmat=nmat,size.thr = .5)})
  # table(Suprathreshold_TF)
  
  #an alternative and faster way:
  Suprathreshold_TF=which(map,arr.ind = TRUE)
  #########
  dd = dist(Suprathreshold_TF)
  hc = hclust(dd, "single")
  # plot(hc)
  # ct = cutree(hc,k=5)
  # pander(table(ct))
  # 
  
  ct = cutree(hc,h=max_dist)
  
  ## sort the cluster names on the basis of their size
  new_cluster_names=rank(table(ct),ties.method = "random")
  ct_new=rep(NA,length(ct))
  for(i in 1:length(new_cluster_names)){
    ct_new[ct==as.numeric(names(new_cluster_names)[i])]=new_cluster_names[i]
  }
  # table(ct,ct_new)
  # table(ct_new)
  ct=ct_new
  rm(ct_new)
  #########
  cluster_map=array(0,dim(map))
  cluster_map[map] =ct
  # pander::pander(table(cluster_map))
  # print(table(cluster_map))
  cluster_map
}
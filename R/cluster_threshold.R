
cluster_threshold <- function(map, max_dist=sqrt(3)){
    ### slower:
    # map=spmT
    # threshold=3.2
    # nmat <- expand.grid(-1:1, -1:1, -1:1)
    # nmat <- nmat[-c(1,3,7,9,14,19,21,25,27), ]
    # system.time(
    # {Suprathreshold_TF = cluster.threshold(spmT>=3.2, nmat=nmat,size.thr = .5)})
    # table(Suprathreshold_TF)
    map = ifelse(is.na(map), FALSE, map)
    #an alternative and faster way:
    Suprathreshold_TF=which(map,arr.ind = TRUE)
    #########
    hc =  hclust.vector(Suprathreshold_TF, "single")
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
    #cluster_map=rep(0,length(map))
    cluster_map[map] =ct
    # pander::pander(table(cluster_map))
    # print(table(cluster_map))
    cluster_map
}
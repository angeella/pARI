
cluster_threshold<- function (map, max_dist = sqrt(3)) 
{
    map = ifelse(is.na(map), FALSE, map)
    Suprathreshold_TF = which(map, arr.ind = TRUE)
    
    dd = dist(Suprathreshold_TF)
    hc = hclust(dd, "single")
    ct = cutree(hc, h = max_dist)
    new_cluster_names = rank(table(ct), ties.method = "random")
    ct_new = rep(NA, length(ct))
    for (i in 1:length(new_cluster_names)) {
        ct_new[ct == as.numeric(names(new_cluster_names)[i])] = new_cluster_names[i]
    }
    ct = ct_new
    rm(ct_new)
    cluster_map = array(0, dim(map))
    cluster_map[map] = ct
    cluster_map
}
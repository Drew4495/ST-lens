# Clustering ----
# |||||||||||||||

## Louvain ----
## ||||||||||||

# Performs Louvain clustering on input low dimensional components.
# - clusternum: The desired number of clusters the user wants to obtain.
# - latent_dat: A d by n matrix of low dimensional components, d is number 
#               of PCs, n is number of spots.
# - knearest: An integers, number of nearest neighbors for KNN graph 
#             construction in louvain clustering.

louvain_clustering = function(clusternum, latent_dat, knearest=100){
  set.seed(1234)
  
  PCvalues = latent_dat
  info.spatial = as.data.frame(t(PCvalues))
  colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))
  knn.norm = FNN::get.knn(as.matrix(t(PCvalues)), k = knearest)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
                                   k=knearest), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = igraph::simplify(nw.norm)
  lc.norm = igraph::cluster_louvain(nw.norm)
  merged <- bluster::mergeCommunities(nw.norm, lc.norm$membership, number=clusternum)
  clusterlabel = as.character(as.integer(as.factor(paste0("cluster",merged))))
  return("cluster_label"=clusterlabel)
}




## Walktrap ----
## |||||||||||||

# This function performs walktrap clustering on input low dimensional components.
# - clusternum: The desired number of clusters the user wants to obtain.
# - latent_dat: A d by n matrix of low dimensional components, d is number of 
#               PCs, n is number of spots.
# - knearest: An integers, number of nearest neighbors for SNN graph 
#             construction in walktrap clustering.

walktrap_clustering = function(clusternum, latent_dat, knearest=100){
  set.seed(1234)
  
  PCvalues = latent_dat
  g <- bluster::makeSNNGraph(as.matrix(t(PCvalues)),k = knearest)
  g_walk <- igraph::cluster_walktrap(g)
  cluster_label_new = as.character(igraph::cut_at(g_walk, no=clusternum))
  
  return("cluster_label"=cluster_label_new)
}



## Refine Cluster 10X ----
## |||||||||||||

# Refines spatial clustering of ST or Visium data.
# - clusterlabels: The cluster label obtained
#   (e.g. from louvain method or walktrap method).
# - location: A n by 2 location matrix of spots.
# - shape: Select shape='hexagon' for Visium data, 'square' for ST data.

refine_cluster_10x = function(clusterlabels, location, shape="square"){
  
  dis_df = as.matrix(dist(location))
  if(shape=="square"){
    num_obs = 4
  }else if(shape == "hexagon"){
    num_obs = 6
  }else{
    print("Select shape='hexagon' for Visium data, 'square' for ST data.")
  }
  refined_pred = clusterlabels
  for(i in 1:length(clusterlabels)){
    nearby_spot_ind = order(dis_df[i,])[1:(num_obs+1)]
    labels_nearby_spot_ind = refined_pred[nearby_spot_ind]  # use updated cluster
    spot_of_interest = refined_pred[i]
    labels_table = table(labels_nearby_spot_ind)
    if( labels_table[spot_of_interest]<num_obs/2 & max(labels_table)>num_obs/2){
      refined_pred[i] = names(labels_table)[which.max(labels_table)]
    }else{
      refined_pred[i] = spot_of_interest
    }
    
  }
  
  return(refined_pred)
}



## fx_CHAOS ----
## |||||||||||||

fx_CHAOS = function(clusterlabel, location){
  # require(parallel)
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  matched_location = scale(matched_location)
  dist_val = rep(0,length(unique(clusterlabel)))
  count = 0
  for(k in unique(clusterlabel)){
    count = count + 1
    location_cluster = matched_location[which(clusterlabel == k),]
    if(length(location_cluster)==2){next}
    #require(parallel)
    results = mclapply(1:dim(location_cluster)[1], fx_1NN, location_in=location_cluster,mc.cores = 5)
    dist_val[count] = sum(unlist(results))
  }
  dist_val = na.omit(dist_val)
  return(sum(dist_val)/length(clusterlabel))
  
}



## fx_1NN ----
## |||||||||||||
#' @import pdist
fx_1NN = function(i,location_in){
  # library(pdist)
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  return(min(line_i))
}



## fx_KNN ----
## |||||||||||||
fx_kNN = function(i,location_in,k,cluster_in){
  #library(pdist)
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  ind = order(line_i)[1:k]
  cluster_use = cluster_in[-i]
  if(sum(cluster_use[ind] != cluster_in[i])>(k/2)){
    return(1)
  }else{
    return(0)
  }
  
}



## fx_PAS ----
## |||||||||||||

fx_PAS = function(clusterlabel, location){
  # require(parallel)
  
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  
  results = mclapply(1:dim(matched_location)[1], fx_kNN, location_in=matched_location,k=10,cluster_in=clusterlabel, mc.cores = 5)
  return(sum(unlist(results))/length(clusterlabel))
}
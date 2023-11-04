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

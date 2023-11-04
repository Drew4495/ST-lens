# Clustering on HR grid ----
# ||||||||||||||||||||||||||

# Assembling data for clustering
data.clustering <- NULL
for(i in 1:nComp){
  loading_HR <- field.eval(grid, loadings_nodes[,i], mesh)
  data.clustering <- cbind(data.clustering, loading_HR)
}

# Clustering
knearest <- round(sqrt(nrow(grid)))
cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
                                         latent_dat = t(data.clustering[,1:nComp_opt]),
                                         knearest = knearest)
if(refine){
  cluster_labels_refined_HR <- refine_cluster_10x(cluster_labels_HR,
                                                  grid,
                                                  shape = refine_type)
} else {
  cluster_labels_refined_HR <- NULL
}

cluster_labels_downsamples <- c()
for(name.l in names.locations){
  
  distances <- dist_point_from_points(locations[name.l,], grid)
  index.closest_point <- which.min(distances)
  if(refine){
    cluster_labels_downsamples[name.l] <- cluster_labels_refined_HR[index.closest_point]
  } else {
    cluster_labels_downsamples[name.l] <- cluster_labels_HR[index.closest_point]
  }
  
}

# Performance index
names.locations_not_unknown <- rownames(locations[true_labels!=7,])
ARI_downsampled <- adjustedRandIndex(true_labels[names.locations_not_unknown,],
                                     cluster_labels_downsamples[names.locations_not_unknown])


## Save clusters ----
## ||||||||||||||||||

save(cluster_labels_HR, cluster_labels_refined_HR, cluster_labels_downsamples,
     ARI_downsampled,
     # Saving options
     file = paste(directory.results, name.dataset, "_clusters_on_HR_grid", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(cluster_labels_HR, cluster_labels_refined_HR, cluster_labels_downsamples,
   ARI_downsampled,
   names.locations_not_unknown)
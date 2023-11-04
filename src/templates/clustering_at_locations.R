# Clustering at locations ----
# ||||||||||||||||||||||||||||

# Assembling data for clustering
data.clustering <- NULL
for(i in 1:nComp){
  data.clustering <- cbind(data.clustering, loadings[,i])
}

# Clustering
cluster_labels <- walktrap_clustering(clusternum = clusternum,
                                      latent_dat = t(data.clustering[,1:nComp_opt]),
                                      knearest = knearest)
if(refine){
  cluster_labels_refined <- refine_cluster_10x(cluster_labels,
                                               locations,
                                               shape = refine_type)
} else {
  cluster_labels_refined <- NULL
}

# Performance index
names.locations_not_unknown <- rownames(locations[true_labels!=7,])

names(cluster_labels) <- names.locations
ARI <- adjustedRandIndex(true_labels[names.locations_not_unknown,],
                         cluster_labels[names.locations_not_unknown])

if(refine){
  names(cluster_labels_refined) <- names.locations
  ARI_refiend <- adjustedRandIndex(true_labels[names.locations_not_unknown,],
                                   cluster_labels_refined[names.locations_not_unknown])
} else {
  ARI_refiend <- NULL
}

## Save clusters ----
## ||||||||||||||||||

save(cluster_labels, cluster_labels_refined,
     ARI, ARI_refiend,
     # Saving options
     file = paste(directory.results, name.dataset, "_clusters_at_locations", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(cluster_labels, cluster_labels_refined,
   ARI, ARI_refiend,
   names.locations_not_unknown)


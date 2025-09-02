library(SpatialPCA)
library(mclust)




get_aligned_labels <- function(clustering_labels, ground_truth) {
  # Initialize an empty data frame to store the ARI scores
  ari_scores <- data.frame(
    ClusterLabel = integer(0),
    GroundTruthLabel = integer(0),
    ARI = numeric(0)
  )
  
  # Calculate ARI for each pairing of clustering label and ground truth label
  for (clust_label in unique(clustering_labels)) {
    for (gt_label in unique(ground_truth)) {
      temp_clust = ifelse(clustering_labels == clust_label, clust_label, 0)
      temp_gt = ifelse(ground_truth == gt_label, gt_label, 0)
      
      ari = adjustedRandIndex(temp_clust, temp_gt)  # Compute ARI
      ari_scores <- rbind(ari_scores, data.frame(ClusterLabel = clust_label, GroundTruthLabel = gt_label, ARI = ari))
    }
  }
  
  # Sort by ARI in descending order
  ari_scores <- ari_scores[order(-ari_scores$ARI),]
  
  # Initialize empty vectors to keep track of final mappings
  final_clust_map <- integer(0)
  final_gt_map <- integer(0)
  
  # Loop through the sorted ARI scores to find optimal unique mappings
  for(i in 1:nrow(ari_scores)) {
    clust_label = ari_scores$ClusterLabel[i]
    gt_label = ari_scores$GroundTruthLabel[i]
    
    if(!(clust_label %in% final_clust_map) && !(gt_label %in% final_gt_map)) {
      final_clust_map <- c(final_clust_map, clust_label)
      final_gt_map <- c(final_gt_map, gt_label)
    }
  }
  
  # Create a mapping from final_clust_map to final_gt_map
  optimal_permutation <- setNames(final_gt_map, final_clust_map)
  
  # Apply this mapping to your original clustering labels
  aligned_clustering_labels <- sapply(clustering_labels, function(label) {
    if(as.character(label) %in% names(optimal_permutation)) {
      return(optimal_permutation[[as.character(label)]])
    } else {
      warning(paste("Label not found, keeping original label:", label))
      return(label)  # Keep the original label
    }
  })
}





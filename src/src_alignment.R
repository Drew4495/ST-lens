#library(SpatialPCA)
library(mclust)




get_aligned_labels <- function(clustering_labels, ground_truth, return_mapping=FALSE) {
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
      
  #Return either aligned labels alone or with mapping
  if (return_mapping) {
    return(list(aligned_clustering_labels = aligned_clustering_labels, 
                optimal_permutation = optimal_permutation))
  } else {
    return(aligned_clustering_labels)
  }
}


align_with_mapping <- function(clustering_labels, mapping){
  aligned_clustering_labels <- sapply(clustering_labels, function(label) {
    if(as.character(label) %in% names(mapping)) {
      return(mapping[[as.character(label)]])
    } else {
      warning(paste("Label not found, keeping original label:", label))
      return(label)  # Keep the original label
    }
  })
}

align_with_nonunique_mapping <- function(clustering_labels, mapping) {
  aligned_clustering_labels <- sapply(clustering_labels, function(label) {
    label <- as.character(label)
    mapped_label <- NULL
    
    for (row in 1:nrow(mapping)) {
      if (label %in% mapping[row, -ncol(mapping)]) {
        mapped_label <- mapping[row, ncol(mapping)]
        break
      }
    }
    
    if (!is.null(mapped_label)) {
      return(mapped_label)
    } else {
      warning(paste("Label not found, keeping original label:", label))
      return(label)  # Keep the original label
    }
  })
  
  return(aligned_clustering_labels)
}









#align_HR_LR(LR_labels, LR_locations, HR_labels, HR_locations, radius, n_spots=100){
#  
#  #initialize variables
#  LR_HR_mapping <- data.frame(LR_labels=unique(LR_labels), HR_labels=rep(0,n))
#  #for each cluster label in LR_locations, sample n spots
#  for (label in LR_HR_mapping$LR_labels){
#    # Get indices of LR locations with the current label
#    label_indices <- which(LR_labels == label)
#    # Decide how many spots to sample based on the number of available labels
#    num_labels <- length(label_indices)
#    num_samples <- if (num_labels >= n_spots) n_spots else num_labels
#    # Sample spots
#    sampled_indices <- sample(label_indices, num_samples, replace = FALSE)
#    sampled_spots <- LR_locations[sampled_indices, ]
#    #for each sampled spot, find all grid nodes in certain radius of sampled spot
#    
#    #return highest frequent grid node label for each sampled spot
#    
#  }
#  
#  
#  #return highest frequent sampled spot for each cluster label
#  #return this mapping
#}



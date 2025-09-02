## ||||||||||||||||
# Alignment Tools ----
## ||||||||||||||||

# Adjusted Rand Index alignment
library(mclust)



## ||||||||||||||||
## Alignment without prior mapping ----
## ||||||||||||||||

get_aligned_labels <- function(clustering_labels, ground_truth, return_mapping=FALSE) {
  #' Align clustering labels to ground truth using Adjusted Rand Index
  #'
  #' This function computes the Adjusted Rand Index (ARI) for each unique pair of clustering
  #' and ground truth labels, determines the optimal alignment based on the highest ARI scores,
  #' and optionally returns the optimal mapping alongside the aligned labels.
  #'
  #' @param clustering_labels Numeric or character vector, labels from clustering.
  #' @param ground_truth Numeric or character vector, ground truth labels.
  #' @param return_mapping Logical, whether to return the mapping as well.
  #' @return A vector of aligned clustering labels or a list containing the labels and mapping.
  #' @import mclust
  
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




## ||||||||||||||||
## Alignment with Mapping ----
## ||||||||||||||||
align_with_mapping <- function(clustering_labels, mapping){
  #' Align clustering labels with a predefined mapping
  #'
  #' This function adjusts clustering labels based on a given mapping of labels,
  #' with warnings for any labels that are not found in the mapping.
  #'
  #' @param clustering_labels Vector, original clustering labels.
  #' @param mapping Named vector, mapping from original to new labels.
  #' @return Vector of aligned clustering labels.
  
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
  #' Align clustering labels with a non-unique mapping table
  #'
  #' Adjusts labels according to a mapping table where each label may be mapped to multiple new labels.
  #' The first applicable mapping encountered is used.
  #'
  #' @param clustering_labels Vector, original clustering labels.
  #' @param mapping Data frame, each row represents a possible mapping for labels.
  #' @return Vector of aligned clustering labels.
  
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


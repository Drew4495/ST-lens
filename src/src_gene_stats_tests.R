
## 2 sample t-test ----
## ||||||||||||||||

between_cluster_ttest <- function(location_labels, labels_of_interest, labeled_gene_by_spots){
  # Load required package
  if (!require("car", character.only = TRUE)) {
    stop("The 'car' package is required but is not installed.")
  }
  
  #initialize variables
  gene.names <- rownames(labeled_gene_by_spots)
  num_genes <- length(gene.names)
  p_values <- numeric(num_genes)
  names(p_values) <- gene.names
  if (is.data.frame(location_labels)){ 
    location_labels <- setNames(unlist(location_labels[, 1, drop = FALSE]), rownames(location_labels))
  }
  
  #subset counts based off labels of interest and non-label
  selected_locations <- names(location_labels[location_labels %in% labels_of_interest])
  selected_counts <- labeled_gene_by_spots[ ,colnames(labeled_gene_by_spots) %in% selected_locations]
  unselected_counts <- labeled_gene_by_spots[ ,!colnames(labeled_gene_by_spots) %in% selected_locations]
  
  for (i in 1:num_genes){
    #Select group1 and group2
    group1 <- selected_counts[i, ]
    group2 <- unselected_counts[i, ]
    
    #Perform Levene Test to decide if variances are equal
    levene_test <- leveneTest(y = c(group1, group2), 
                              group = factor(rep(1:2, c(length(group1), length(group2)))))
    var.equal <- levene_test$`Pr(>F)`[1] > 0.05
    
    # Check if levene_test$p.value exists
    if (!is.null(levene_test$`Pr(>F)`[1])) {
      var.equal <- levene_test$`Pr(>F)`[1] > 0.05
    } else {
      var.equal <- FALSE
      warning(paste("Levene's test p-value not available for gene", gene.names[i]))
    }
    
   
    #Perform t-test
    t_test_result <- t.test(group1, group2, var.equal = var.equal)
    p_values[i] <- t_test_result$p.value
  }
  
  #Apply FDR correction
  p_adj <- p.adjust(p_values, method = 'BH')
  return(p_adj)
}



## ANOVA ----
## ||||||||||||||||



## get_gene_means_by_cluster ----
## ||||||||||||||||

get_gene_means_by_cluster <- function(location_labels, labeled_gene_by_spots, genes_of_interest){
  
  #initialize variables
  labels <- sort(unique(location_labels), decreasing=FALSE)
  labeled_gene_by_spots <- labeled_gene_by_spots[genes_of_interest,]
  mean_by_cluster_df <- data.frame(matrix(NA, nrow=dim(labeled_gene_by_spots)[1], ncol=length(labels)))
  colnames(mean_by_cluster_df) <- labels
  rownames(mean_by_cluster_df) <- genes_of_interest
  
  #loop through locations by label and calculate mean
  for (lab in labels){
    selected_locations <- names(location_labels[location_labels %in% lab])
    selected_counts <- labeled_gene_by_spots[ ,colnames(labeled_gene_by_spots) %in% selected_locations]
      
    row_means <- apply(selected_counts, 1, mean)
    mean_by_cluster_df[ ,lab] <- row_means
    
  }
  
  return(mean_by_cluster_df)
}




# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Clustering SlideseqCerebellum dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()





## ||||||||||||||||||||
# libraries and functions ----
## ||||||||||||||||||||

## plots
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(patchwork)

## Spatial data
library(sp)
library(sf)
library(rmapshaper)
library(spatstat)

## fdaPDE
library(fdaPDE)
library(fdaPDE2)

## ARI performance index
library(mclust)

## parallel
library(foreach)
library(doParallel)

## Other utils
library(glue)

## functions
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/plots.R")
source("src/clustering.R")
source("src/src_alignment.R")



## ||||||||||||||||||||
# Set globals ----
## ||||||||||||||||||||

## dataset name
name_dataset <- "SlideseqCerebellum"

cat.script_title(glue("Clustering {name_dataset} dataset"))

# global variables 
cat.section_title("Global variables")


## directories
path_data <- paste("data/", name_dataset, "/", sep = "")
path_images <- paste("images/", name_dataset, "/", sep = "")
mkdir(c(path_images))
path_images <- paste("images/", name_dataset, "/clustering/", sep = "")
mkdir(c(path_images))
path_results <- paste("results/", name_dataset, "/", sep = "")
mkdir(c(path_results))

## figures dimensions
figure_width <- 8
figure_height <- 8




#=============================================================================#



## ||||||||||||||||||||
# Load in Data ----
## ||||||||||||||||||||

load(paste(path_data, "mesh.RData", sep = ""))
load(paste(path_data, "analyzed_data.RData", sep = ""))
load(paste(path_results, "fPCA.RData", sep = ""))
load(paste(path_results, "nComp_opt.RData", sep = ""))

rm(loadings_locs_NotNormalized, loadings_HR_NotNormalized)




#=============================================================================#



## ||||||||||||||||||||
# Define functions for aligning labels to ground truth----
## ||||||||||||||||||||

align_all_list_labels <- function(old_label_list, gt_labels) {
  new_label_list <- list()  # Initialize the new list within the function
  for (label_name in names(old_label_list)) {
    label_vec <- old_label_list[[label_name]]
    aligned_labels <- get_aligned_labels(clustering_labels = label_vec, ground_truth = gt_labels)
    new_label_list[[label_name]] <- aligned_labels
  }
  return(new_label_list)  # Return the updated list
}


replot_and_aggregate_labels <- function(cluster_labels_HR_list, cluster_labels_list, grid, 
                                        locations, dir_output, figure_width=8, figure_height=8, 
                                        do_aggregate = TRUE, LR_size = 2) {
  dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)
  
  HR_plots <- list()
  LR_plots <- list()
  knn_values <- c()
  
  for (key in names(cluster_labels_HR_list)) {
    knearest <- as.numeric(gsub("knearest_(\\d+).*", "\\1", key))
    knn_values <- c(knn_values, knearest)
    
    # Plot HR aligned labels
    plot_HR <- plot.field_points(grid, cluster_labels_HR_list[[key]], 
                                 colormap = "H", size = 0.5, discrete = TRUE) +
      ggtitle(paste("HR Clustering - knn:", knearest)) + 
      standard_plot_settings_fields()
    ggsave(paste(dir_output, "HR_clustering_aligned_k", knearest , ".pdf", sep = ""),
           plot = plot_HR, width = figure_width, height = figure_height, dpi = "print", units = "cm")
    
    # Store HR plot for aggregation
    HR_plots[[key]] <- plot_HR + ggtitle(paste(knearest)) +
      theme(legend.position = "none")
    
    # Plot original location aligned labels
    plot_LR <- plot.field_points(locations, cluster_labels_list[[key]], 
                                 colormap = "H", size = LR_size, discrete = TRUE) +
      theme(legend.position = "none") +
      ggtitle(paste("Clustering at locations - knn:", knearest)) + 
      standard_plot_settings_fields()
    ggsave(paste(dir_output, "clustering_aligned_k", knearest , ".pdf", sep = ""),
           plot = plot_LR, width = figure_width, height = figure_height, dpi = "print", units = "cm")
    
    # Store LR plot for aggregation
    LR_plots[[key]] <- plot_LR + ggtitle(paste(knearest)) + 
      theme(legend.position = "none")
  }
  
  # Perform aggregation if required
  if (do_aggregate) {
    # Combine HR plots in a 5x6 grid
    combined_HR_plot <- wrap_plots(HR_plots, ncol = 6, nrow = 5)
    ggsave(paste(dir_output, "aggregated_KNN_HR_plot.pdf", sep = ""),
           plot = combined_HR_plot, width = figure_width * 6, height = figure_height * 5, dpi = "print", units = "cm")
    
    # Combine LR plots in a 5x6 grid
    combined_LR_plot <- wrap_plots(LR_plots, ncol = 6, nrow = 5)
    ggsave(paste(dir_output, "aggregated_KNN_LR_plot.pdf", sep = ""),
           plot = combined_LR_plot, width = figure_width * 6, height = figure_height * 5, dpi = "print", units = "cm")
  }
}





#=============================================================================#



## ||||||||||||||||||||
# clustering ----
## ||||||||||||||||||||


## Prep for clustering ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## clustering grid
grid_step <- 25
grid <- square_grid(SpatialPoints(locations)@bbox, grid_step, seed_point = seed_point)
grid <- grid[!is.na(over(SpatialPoints(grid), lattice$domain)),]

## grid plot
plot <- ggplot() +
  standard_plot_settings_fields() + ggtitle("HR grid") +
  geom_sf(data = st_as_sf(SpatialPoints(grid)), color = "black", size = 0.25) +
  geom_sf(data = st_as_sf(SpatialPoints(locations)), color = "red", size = 0.25)
ggsave(paste(path_images, "1_HR_grid.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 1. nsc = 3 ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Define variables
#loadings <- loadings
#mean <- mean

## hyper-parameters
nsc <- 3
clusternum <- 9 
#knearest_vect <- seq(10, 300, by = 10)  # 200 # round(sqrt(nrow(grid)))
knearest <- round(sqrt(nrow(grid)))
knearest_vect <- c(knearest)

## Make clustering directory
dir_clustering_images <- file.path(path_images, glue("GCV_NoMean/nonaligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images, recursive = TRUE, showWarnings = FALSE)
dir_clustering_images_aligned <- file.path(path_images, "GCV_NoMean/", glue("aligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images_aligned, recursive = TRUE, showWarnings = FALSE)

## mean at the clustering grid
mean_clustering <- evaluate_field(grid, mean, mesh)

## functional components evaluation at the clustering grid
loadings_clustering <- NULL
for(i in 1:n_comp) {
  loadings_clustering <- cbind(loadings_clustering, evaluate_field(grid, loadings[, i], mesh))
}

## Set this variable to TRUE if you want to generate and save plots
is_plot <- TRUE

# Initialize lists to store the cluster labels and ARI values
cluster_labels_HR_list <- list()
cluster_labels_list <- list()
#ARI_list <- list()

for (knearest in knearest_vect) {
  
  ## data for clustering
  plus_mean <- 0
  label_mean <- ""
  data_clustering <- loadings_clustering
  include_mean <- FALSE
  if(include_mean){
    data_clustering <- cbind(mean_clustering, loadings_clustering)
    plus_mean <- 1
    label_mean <- "_mean"
  }
  
  ## clustering
  #cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
  #                                         latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
  #                                         knearest = knearest)
  cluster_labels_HR <- louvain_clustering(clusternum = clusternum,
                                          latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
                                          knearest = knearest)
  
  if (is_plot) {
    ## plot clustering HR
    plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 0.2, discrete = TRUE) +
      ggtitle("HR Clustering") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "HR_clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")
  }
  
  ## computing clustering at the original locations
  cluster_labels <- c()
  for(names_location in names_locations){
    distances <- dist_point_from_points(locations[names_location,], grid)
    index.closest_point <- which.min(distances)
    cluster_labels[names_location] <- cluster_labels_HR[index.closest_point]
  }
  
  if (is_plot) {
    ## plot clustering at locations
    plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 0.1, discrete = TRUE) +
      ggtitle("Clustering at locations") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
  }
  
  #Store results
  key <- paste0("knearest_", knearest, label_mean)
  cluster_labels_HR_list[[key]] <- cluster_labels_HR
  cluster_labels_list[[key]] <- cluster_labels
  #ARI_list[[key]] <- ARI
}


#### Save variables ----
cluster_labels_HR_list_GCV_NoMean <- cluster_labels_HR_list
cluster_labels_list_GCV_NoMean <- cluster_labels_list
#ARI_list_GCV_NoMean <- ARI_list
filepath <- paste(path_results, glue("fPCA_clustering_results_GCV_NoMean_{nsc}.RData"), sep = "")
save(cluster_labels_HR_list_GCV_NoMean, cluster_labels_list_GCV_NoMean,
     file = filepath
)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 2. nsc = 4 ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Define variables
#loadings <- loadings
#mean <- mean

## hyper-parameters
nsc <- 4
clusternum <- 9 
#knearest_vect <- seq(10, 300, by = 10)  # 200 # round(sqrt(nrow(grid)))
knearest <- round(sqrt(nrow(grid)))
knearest_vect <- c(knearest)

## Make clustering directory
dir_clustering_images <- file.path(path_images, glue("GCV_NoMean/nonaligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images, recursive = TRUE, showWarnings = FALSE)
dir_clustering_images_aligned <- file.path(path_images, "GCV_NoMean/", glue("aligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images_aligned, recursive = TRUE, showWarnings = FALSE)

## mean at the clustering grid
mean_clustering <- evaluate_field(grid, mean, mesh)

## functional components evaluation at the clustering grid
loadings_clustering <- NULL
for(i in 1:n_comp) {
  loadings_clustering <- cbind(loadings_clustering, evaluate_field(grid, loadings[, i], mesh))
}

## Set this variable to TRUE if you want to generate and save plots
is_plot <- TRUE

# Initialize lists to store the cluster labels and ARI values
cluster_labels_HR_list <- list()
cluster_labels_list <- list()
#ARI_list <- list()

for (knearest in knearest_vect) {
  
  ## data for clustering
  plus_mean <- 0
  label_mean <- ""
  data_clustering <- loadings_clustering
  include_mean <- FALSE
  if(include_mean){
    data_clustering <- cbind(mean_clustering, loadings_clustering)
    plus_mean <- 1
    label_mean <- "_mean"
  }
  
  ## clustering
  #cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
  #                                         latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
  #                                         knearest = knearest)
  cluster_labels_HR <- louvain_clustering(clusternum = clusternum,
                                          latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
                                          knearest = knearest)
  
  if (is_plot) {
    ## plot clustering HR
    plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 0.2, discrete = TRUE) +
      ggtitle("HR Clustering") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "HR_clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")
  }
  
  ## computing clustering at the original locations
  cluster_labels <- c()
  for(names_location in names_locations){
    distances <- dist_point_from_points(locations[names_location,], grid)
    index.closest_point <- which.min(distances)
    cluster_labels[names_location] <- cluster_labels_HR[index.closest_point]
  }
  
  if (is_plot) {
    ## plot clustering at locations
    plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 0.1, discrete = TRUE) +
      ggtitle("Clustering at locations") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
  }
  
  #Store results
  key <- paste0("knearest_", knearest, label_mean)
  cluster_labels_HR_list[[key]] <- cluster_labels_HR
  cluster_labels_list[[key]] <- cluster_labels
  #ARI_list[[key]] <- ARI
}


#### Save variables ----
cluster_labels_HR_list_GCV_NoMean <- cluster_labels_HR_list
cluster_labels_list_GCV_NoMean <- cluster_labels_list
#ARI_list_GCV_NoMean <- ARI_list
filepath <- paste(path_results, glue("fPCA_clustering_results_GCV_NoMean_{nsc}.RData"), sep = "")
save(cluster_labels_HR_list_GCV_NoMean, cluster_labels_list_GCV_NoMean,
     file = filepath
)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 2. nsc = 5 (nComp_opt) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Define variables
#loadings <- loadings
#mean <- mean

## hyper-parameters
nsc <- 5
clusternum <- 9 
#knearest_vect <- seq(10, 300, by = 10)  # 200 # round(sqrt(nrow(grid)))
knearest <- round(sqrt(nrow(grid)))
knearest_vect <- c(knearest)

## Make clustering directory
dir_clustering_images <- file.path(path_images, glue("GCV_NoMean/nonaligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images, recursive = TRUE, showWarnings = FALSE)
dir_clustering_images_aligned <- file.path(path_images, "GCV_NoMean/", glue("aligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images_aligned, recursive = TRUE, showWarnings = FALSE)

## mean at the clustering grid
mean_clustering <- evaluate_field(grid, mean, mesh)

## functional components evaluation at the clustering grid
loadings_clustering <- NULL
for(i in 1:n_comp) {
  loadings_clustering <- cbind(loadings_clustering, evaluate_field(grid, loadings[, i], mesh))
}

## Set this variable to TRUE if you want to generate and save plots
is_plot <- TRUE

# Initialize lists to store the cluster labels and ARI values
cluster_labels_HR_list <- list()
cluster_labels_list <- list()
#ARI_list <- list()

for (knearest in knearest_vect) {
  
  ## data for clustering
  plus_mean <- 0
  label_mean <- ""
  data_clustering <- loadings_clustering
  include_mean <- FALSE
  if(include_mean){
    data_clustering <- cbind(mean_clustering, loadings_clustering)
    plus_mean <- 1
    label_mean <- "_mean"
  }
  
  ## clustering
  #cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
  #                                         latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
  #                                         knearest = knearest)
  cluster_labels_HR <- louvain_clustering(clusternum = clusternum,
                                          latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
                                          knearest = knearest)
  
  if (is_plot) {
    ## plot clustering HR
    plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 0.2, discrete = TRUE) +
      ggtitle("HR Clustering") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "HR_clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")
  }
  
  ## computing clustering at the original locations
  cluster_labels <- c()
  for(names_location in names_locations){
    distances <- dist_point_from_points(locations[names_location,], grid)
    index.closest_point <- which.min(distances)
    cluster_labels[names_location] <- cluster_labels_HR[index.closest_point]
  }
  
  if (is_plot) {
    ## plot clustering at locations
    plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 0.1, discrete = TRUE) +
      ggtitle("Clustering at locations") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
  }
  
  #Store results
  key <- paste0("knearest_", knearest, label_mean)
  cluster_labels_HR_list[[key]] <- cluster_labels_HR
  cluster_labels_list[[key]] <- cluster_labels
  #ARI_list[[key]] <- ARI
}


#### Save variables ----
cluster_labels_HR_list_GCV_NoMean <- cluster_labels_HR_list
cluster_labels_list_GCV_NoMean <- cluster_labels_list
#ARI_list_GCV_NoMean <- ARI_list
filepath <- paste(path_results, glue("fPCA_clustering_results_GCV_NoMean_{nsc}.RData"), sep = "")
save(cluster_labels_HR_list_GCV_NoMean, cluster_labels_list_GCV_NoMean,
     file = filepath
)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 2. nsc = 6 (nComp_opt) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Define variables
#loadings <- loadings
#mean <- mean

## hyper-parameters
nsc <- 6
clusternum <- 9 
#knearest_vect <- seq(10, 300, by = 10)  # 200 # round(sqrt(nrow(grid)))
knearest <- round(sqrt(nrow(grid)))
knearest_vect <- c(knearest)

## Make clustering directory
dir_clustering_images <- file.path(path_images, glue("GCV_NoMean/nonaligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images, recursive = TRUE, showWarnings = FALSE)
dir_clustering_images_aligned <- file.path(path_images, "GCV_NoMean/", glue("aligned_labels/NSC_{nsc}/"))
dir.create(dir_clustering_images_aligned, recursive = TRUE, showWarnings = FALSE)

## mean at the clustering grid
mean_clustering <- evaluate_field(grid, mean, mesh)

## functional components evaluation at the clustering grid
loadings_clustering <- NULL
for(i in 1:n_comp) {
  loadings_clustering <- cbind(loadings_clustering, evaluate_field(grid, loadings[, i], mesh))
}

## Set this variable to TRUE if you want to generate and save plots
is_plot <- TRUE

# Initialize lists to store the cluster labels and ARI values
cluster_labels_HR_list <- list()
cluster_labels_list <- list()
#ARI_list <- list()

for (knearest in knearest_vect) {
  
  ## data for clustering
  plus_mean <- 0
  label_mean <- ""
  data_clustering <- loadings_clustering
  include_mean <- FALSE
  if(include_mean){
    data_clustering <- cbind(mean_clustering, loadings_clustering)
    plus_mean <- 1
    label_mean <- "_mean"
  }
  
  ## clustering
  #cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
  #                                         latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
  #                                         knearest = knearest)
  cluster_labels_HR <- louvain_clustering(clusternum = clusternum,
                                          latent_dat = t(data_clustering[,1:(plus_mean + nsc)]),
                                          knearest = knearest)
  
  if (is_plot) {
    ## plot clustering HR
    plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 0.2, discrete = TRUE) +
      ggtitle("HR Clustering") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "HR_clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")
  }
  
  ## computing clustering at the original locations
  cluster_labels <- c()
  for(names_location in names_locations){
    distances <- dist_point_from_points(locations[names_location,], grid)
    index.closest_point <- which.min(distances)
    cluster_labels[names_location] <- cluster_labels_HR[index.closest_point]
  }
  
  if (is_plot) {
    ## plot clustering at locations
    plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 0.1, discrete = TRUE) +
      ggtitle("Clustering at locations") + standard_plot_settings_fields() + theme(legend.position="none")
    ggsave(paste(dir_clustering_images, "clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
           plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
  }
  
  #Store results
  key <- paste0("knearest_", knearest, label_mean)
  cluster_labels_HR_list[[key]] <- cluster_labels_HR
  cluster_labels_list[[key]] <- cluster_labels
  #ARI_list[[key]] <- ARI
}


#### Save variables ----
cluster_labels_HR_list_GCV_NoMean <- cluster_labels_HR_list
cluster_labels_list_GCV_NoMean <- cluster_labels_list
#ARI_list_GCV_NoMean <- ARI_list
filepath <- paste(path_results, glue("fPCA_clustering_results_GCV_NoMean_{nsc}.RData"), sep = "")
save(cluster_labels_HR_list_GCV_NoMean, cluster_labels_list_GCV_NoMean,
     file = filepath
)

#=============================================================================#






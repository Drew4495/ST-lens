# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Analysis HER2 dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

# libraries ----

## plots
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)

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


# functions ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/plots.R")
source("src/clustering.R")


##############################################################################


cat.script_title("Clustering HER2 dataset")


# global variables ----

cat.section_title("Global variables")

## dataset name
name_dataset <- "HER2"

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


# load data ----
load(paste(path_data, "mesh.RData", sep = ""))
load(paste(path_data, "analyzed_data.RData", sep = ""))
load(paste(path_results, "fPCA.RData", sep = ""))


# clustering ----

## clustering grid
grid_step <- 1/3
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))
grid <- square_grid(SpatialPoints(locations)@bbox, grid_step, seed_point = seed_point)
grid <- grid[!is.na(over(SpatialPoints(grid), lattice$domain)),]

## grid plot
plot <- ggplot() +
  standard_plot_settings_fields() + ggtitle("HR grid") +
  geom_sf(data = st_as_sf(SpatialPoints(grid)), color = "black", size = 0.25) +
  geom_sf(data = st_as_sf(SpatialPoints(locations)), color = "red", size = 0.25)
ggsave(paste(path_images, "1_HR_grid.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

# ## mean at the clustering grid
# mean_clustering <- evaluate_field(grid, mean, mesh)

## functional components evaluation at the clustering grid
loadings_clustering <- NULL
for(i in 1:n_comp) {
  loadings_clustering <- cbind(loadings_clustering, evaluate_field(grid, loadings[, i], mesh))
}

## hyper-parameters
nComp_opt <- 3
clusternum <- 7
knearest_vect <- seq(10, 300, by = 10) # 200 # round(sqrt(nrow(grid)))

## clustering (in parallel)
parallel::detectCores()
n.cores <- parallel::detectCores() - 4
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
foreach(
  knearest = knearest_vect,
  .combine = 'c',
  .packages = c("ggplot2", "viridis", "sf", "sp", "mclust")
) %dopar% {

  ## data for clustering
  plus_mean <- 0
  label_mean <- ""
  data_clustering <- loadings_clustering
  if(FALSE){
    data_clustering <- cbind(mean_clustering, loadings_clustering)
    plus_mean <- 1
    label_mean <- "_mean"
  }
  
  ## clustering
  cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
                                           latent_dat = t(data_clustering[,1:(nComp_opt + plus_mean)]),
                                           knearest = knearest)
  
  ## plot clustering HR
  plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 0.5, discrete = TRUE) +
    ggtitle("HR Clustering") + standard_plot_settings_fields()
  ggsave(paste(path_images, "HR_clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
         plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
  
  
  ## computing clustering at the original locations
  cluster_labels <- c()
  for(names_location in names_locations){
    distances <- dist_point_from_points(locations[names_location,], grid)
    index.closest_point <- which.min(distances)
    cluster_labels[names_location] <- cluster_labels_HR[index.closest_point]
  }
  
  ## plot clustering at locations
  plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 2, discrete = TRUE) +
    ggtitle("Clustering at locations") + standard_plot_settings_fields()
  ggsave(paste(path_images, "clustering", label_mean, "_k", knearest , ".pdf", sep = ""),
         plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")
  
  # Performance index
  ARI <- adjustedRandIndex(true_labels$true_label, cluster_labels)
  
  ARI

}



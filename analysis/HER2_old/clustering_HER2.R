# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Analysis HER2 dataset - old version %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


cat.script_title("Clustering HER2 dataset - old version")


# global variables ----

cat.section_title("Global variables")

## dataset name
name_dataset <- "HER2"

## directories
path_data <- paste("data/", name_dataset, "/", sep = "")
path_images <- paste("images/", name_dataset, "_old/", sep = "")
mkdir(c(path_images))
path_images <- paste("images/", name_dataset, "_old/clustering/", sep = "")
mkdir(c(path_images))
path_results <- paste("results/", name_dataset, "/", sep = "")
mkdir(c(path_results))

## figures dimensions
figure_width <- 8
figure_height <- 8


# load data ----
load(paste(path_data, "mesh_old.RData", sep = ""))
load(paste(path_data, "analyzed_data_old.RData", sep = ""))
load(paste(path_results, "fPCA_old.RData", sep = ""))


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

## functional components evaluation at the clustering grid
loadings_clustering <- NULL
for(i in 1:n_comp) {
  loadings_clustering <- cbind(loadings_clustering, evaluate_field(grid, loadings[, i], mesh))
}

## plot components tile
plot_list <- list()
limits = range(loadings)
for(i in 1:n_comp) {
  plot_list[[i]] <- plot.field_tile(grid, loadings_clustering[, i], limits = limits) +
    standard_plot_settings_fields() +
    ggtitle(paste("w", i, sep = ""))
}
plot <- arrangeGrob(grobs = plot_list, nrow = ceiling(sqrt(n_comp)))
grid.arrange(plot)
ggsave(paste(path_images, "2_functional_components.pdf", sep = ""),
       plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")


## hyper-parameters
nComp_opt <- 3
clusternum <- 7
knearest <- 200

## data for clustering
data_clustering <- loadings_clustering

## clustering
cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
                                         latent_dat = t(data_clustering[,1:(nComp_opt)]),
                                         knearest = knearest)

## plot clustering HR
plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 0.5, discrete = TRUE) +
  ggtitle("HR Clustering") + standard_plot_settings_fields()
ggsave(paste(path_images, "3_HR_clustering", "_k", knearest , ".pdf", sep = ""),
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
ggsave(paste(path_images, "4_clustering", "_k", knearest , ".pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## performance index

## As we were doing
names.locations_not_unknown <- rownames(locations[true_labels!=7,])
ARI_7 <- adjustedRandIndex(true_labels[names.locations_not_unknown,],
                           cluster_labels[names.locations_not_unknown])
ARI_7

## As SpatialPCA is doing
ARI <- adjustedRandIndex(true_labels$true_label, cluster_labels)
ARI


# saving results ----

## clustering results
save(
  ARI,
  grid,
  cluster_labels_HR, cluster_labels,
  file = paste(path_results, "fPCA_old.RData", sep = "")
)

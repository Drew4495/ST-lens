# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Analysis HER2 dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
PLOT <- TRUE
fdaPDE <- FALSE

# setwd("~/ST-lens")
source("src/cat_utilities.R")

cat.script_title("Analysis HER2 dataset")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

source("src/templates/libraries.R")


# ||||||||||||||
# Functions ----
# ||||||||||||||

cat.section_title("Functions")

source("src/plot_utilities.R")
source("src/plots.R")

source("src/meshing.R")
source("src/fdaPDE.R")
source("src/clustering.R")


# |||||||||||||||||||||
# Global variables ----
# |||||||||||||||||||||

cat.section_title("Global variables")

# Directories
directory.initial_data <- "data/HER2/"
directory.results <- "results/HER2/"
directory.images <- "images/HER2/"
name.dataset <- "HER2"

# Code flow control
RUN <- list()
RUN[["Pre-Processing"]] <- FALSE
RUN[["Mesh Generation"]] <- FALSE
RUN[["Mean estimation"]] <- FALSE
RUN[["Optimal nComp selection"]] <- FALSE
RUN[["Optimal components"]] <- FALSE
RUN[["Clustering at locations"]] <- FALSE
RUN[["Clustering on HR grid"]] <- TRUE


# |||||||||||||||||||
# Initialization ----
# |||||||||||||||||||

cat.section_title("Initialization")

source("src/templates/initialization.R")


# |||||||||
# Data ----
# |||||||||

cat.section_title("Data")

load(paste(directory.initial_data, name.dataset, ".RData", sep = ""))

# Data content
# - locations:      (data.frame)    #locations x 2
# - counts:         (dgCMatrix)     #genes x #locations
# - true_labels:    (data.frame)    #locations x 1  (NULL if not available)

# Save initial data
locations.initial <- locations
counts.initial <- counts
true_labels.initial <- true_labels
names.locations.initial <- rownames(locations)
names.genes.initial <- rownames(counts)


# |||||||||||||||||||
# Pre-Processing ----
# |||||||||||||||||||

sparkversion <- "spark" # "sparkX"
numCores_spark <- 1
number.genes <- 3000
min.loctions <- 20
min.features <- 20

# Pre-processing
tic()
if(RUN[["Pre-Processing"]])
  source("src/templates/pre_processing.R")
elapsed <- toc(quiet = TRUE)
times[["Pre-Processing"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load processed data
load(paste(directory.results, name.dataset, "_processed.RData", sep = ""))

# Processed data
locations <- locations.significant
counts <- counts.significant
counts.normalized <- counts.normalized
names.locations <- rownames(locations)
names.genes <- rownames(counts)

# Stats
if(VERBOSE){
  cat("\nStats")
  cat(paste("\n- Initial number of locations:", nrow(locations.initial)))
  cat(paste("\n- Final number of locations:", nrow(locations)))
  cat(paste("\n- Initial number of genes: ", nrow(counts.initial)))
  cat(paste("\n- Final number of genes: ", nrow(counts)))
}

# Clean
rm(locations.significant, counts.significant)


# ||||||||||||||||||||
# Mesh Generation ----
# ||||||||||||||||||||

cat.section_title("Mesh Generation")

# Data exploration
if(PLOT){
  ggplot() +
    standard_plot_settings() + 
    xlab("x") + ylab("y") + ggtitle("Initial Locations") +
    geom_sf(data = st_as_sf(SpatialPoints(locations)), color = "black", size = 1)
}

# Hyper-parameters
h <- 1
bbox <- NULL
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))
type <- "square"
simplification <- 0.15
remove_holes <- FALSE
minimum_area_hole <- NULL
simplification_hole <- 1
maximum_area <- 0.3

# Mesh generation
tic()
if(RUN[["Mesh Generation"]])
  source("src/templates/mesh_generation.R")
elapsed <- toc(quiet = TRUE)
times[["Mesh Generation"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load generated data
load(paste(directory.results, name.dataset, "_mesh.RData", sep = ""))

# Processed data
locations <- locations.final
names.locations <- rownames(locations)
mesh <- mesh
lattice <- lattice

# Stats
if(VERBOSE){
  cat("\nStats")
  cat(paste("\n- Number of locations:", nrow(locations)))
  cat(paste("\n- Number of elements: ", nrow(mesh$triangles)))
  cat(paste("\n- Number of nodes: ", nrow(mesh$nodes)))
}

# Plot final locations
if(PLOT){
  plot <- plot.final_locations(SpatialPoints(locations.initial),
                               SpatialPoints(locations), lattice)
  plot <- plot + xlab("") + ylab("") + ggtitle("Final locations")
  ggsave(paste(directory.images, "final_locations.jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}

# Plot mesh
if(PLOT){
  plot <- plot.fdaPDE_mesh(mesh)
  plot <- plot + xlab("") + ylab("") + ggtitle("Mesh")
  ggsave(paste(directory.images, "mesh.jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}

# Clean
rm(locations.final)


# |||||||||||||||||||||||
# Data decomposition ----
# |||||||||||||||||||||||

cat.section_title("Data Decomposition")

# Update counts about eventually discarded locations:
counts <- counts[, names.locations]
counts.normalized <- counts.normalized[, names.locations]
true_labels <- data.frame(true_label = true_labels.initial[names.locations, ])
row.names(true_labels) <- names.locations

# HR grid
grid <- square_grid(SpatialPoints(locations)@bbox, 1/3, seed_point = seed_point)
grid <- grid[!is.na(over(SpatialPoints(grid), lattice$domain)),]

# Plot
if(PLOT){
  plot <- ggplot() +
    standard_plot_settings() + 
    xlab("x") + ylab("y") + ggtitle("HR grid") +
    geom_sf(data = st_as_sf(SpatialPoints(grid)), color = "black", size = 0.25) +
    geom_sf(data = st_as_sf(SpatialPoints(locations)), color = "red", size = 0.25)
  ggsave(paste(directory.images, "HR_grid.jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}


## Mean estimation ----
## ||||||||||||||||||||

cat.subsection_title("Mean estimation")

# Hyper-parameters
lambdas <- 10^seq(-9, 2, by = 1)

# Mean estimation
tic()
if(RUN[["Mean estimation"]])
  source("src/templates/mean_estimation.R")
elapsed <- toc(quiet = TRUE)
times[["Mean estimation"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load generated data
load(paste(directory.results, name.dataset, "_mean.RData", sep = ""))

# Processed data
counts.mean <- counts.mean
counts.mean_nodes <- counts.mean_nodes
counts.centered <- counts.centered
lambda_opt <- lambda_opt

# Plot
if(PLOT){
  counts.mean_HR <- field.eval(grid, counts.mean_nodes, mesh)
  plot <- plot.field_tile(grid, counts.mean_HR, colormap = "D") +
    ggtitle("Mean genes expression") + xlab("") + ylab("")
  ggsave(paste(directory.images, "mean.jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}


## Optimal nComp selection ----
## ||||||||||||||||||||||||||||

# Hyper-parameters
lambda <- lambda_opt
nComp <- 20
nComp_opt <- 20 # 3 (20 because we want to compute ARI for all the components)

tic()
if(RUN[["Optimal nComp selection"]])
  source("src/templates/nComp_selection.R")
elapsed <- toc(quiet = TRUE)
times[["Optimal nComp selection"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load generated data
load(paste(directory.results, name.dataset, "_nComp_selection.RData", sep = ""))

# Processed data
nComp_opt <- nComp_opt
residuals_norm <- residuals_norm
scores <- scores
loadings <- loadings
loadings_nodes <- loadings_nodes

# Plot components HR all
if(PLOT){
  loadings_HR <- NULL
  for(h in 1:nComp){
    loadings_HR <- cbind(loadings_HR, field.eval(grid, loadings_nodes[,h], mesh))
  }
  plot <- plot.components(grid, loadings_HR, type = "tile")
  ggsave(paste(directory.images, "components_all_HR.jpg", sep = ""),
         plot = plot, width = 3*5, height = 3*4, dpi = 200)
}

# Plot nComp selection
if(PLOT){
  plot <- plot.nComp_selection(residuals_norm, nComp, nComp_opt)
  ggsave(paste(directory.images, "nComp_selection.jpg", sep = ""),
         plot = plot, width = 9, height = 6, dpi = 200)
}

# Clean
rm(residuals_norm, scores, loadings, loadings_nodes)


## Optimal components ----
## |||||||||||||||||||||||

# Hyper-parameters
lambdas <- 10^seq(-9, 2, by = 0.2)
nComp <- nComp_opt

# Optimal components
tic()
if(RUN[["Optimal components"]])
  source("src/templates/optimal_components.R")
elapsed <- toc(quiet = TRUE)
times[["Optimal components"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load generated data
load(paste(directory.results, name.dataset, "_components.RData", sep = ""))

# Processed data
scores <- scores
loadings <- loadings
loadings_nodes <- loadings_nodes

# Plot components HR
if(PLOT){
  loadings_HR <- NULL
  for(h in 1:nComp){
    loadings_HR <- cbind(loadings_HR, field.eval(grid, loadings_nodes[,h], mesh))
  }
  plot <- plot.components(grid, loadings_HR, type = "tile")
  ggsave(paste(directory.images, "components_selected_HR.jpg", sep = ""),
         plot = plot, width = 3*5, height = 3*4, dpi = 200)
}


# Clustering ----
# |||||||||||||||

# Hyper-parameters
nComp_opt <- 3 
clusternum <- 6
knearest <- 200 # round(sqrt(nrow(grid)))

tic()

if(RUN[["Clustering on HR grid"]]){
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  # Assembling data for clustering
  data.clustering <- NULL
  for(h in 1:max(nComp_opt)){
    loading_HR <- field.eval(grid, loadings_nodes[,h], mesh)
    data.clustering <- cbind(data.clustering, loading_HR)
  }
  
  # Clustering
  cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
                                           latent_dat = t(data.clustering[,1:nComp_opt]),
                                           knearest = knearest)
  
  cluster_labels <- c()
  for(name.l in names.locations){
    
    distances <- dist_point_from_points(locations[name.l,], grid)
    index.closest_point <- which.min(distances)
    cluster_labels[name.l] <- cluster_labels_HR[index.closest_point]
    
  }
  
  # Performance index
  names.locations_not_unknown <- rownames(locations[true_labels!=7,])
  ARI <- adjustedRandIndex(true_labels[names.locations_not_unknown,],
                           cluster_labels[names.locations_not_unknown])
  
  # Save clusters
  save(cluster_labels_HR, cluster_labels,
       ARI, nComp_opt, grid,
       # Saving options
       file = paste(directory.results, name.dataset, "_clusters_on_HR_grid", "_nComp", nComp_opt, ".RData", sep = ""))
  
  # Clean
  rm(cluster_labels_HR, cluster_labels,
     ARI,
     names.locations_not_unknown)
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
}

elapsed <- toc(quiet = TRUE)
times[["Clustering on HR grid"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load generated data
load(paste(directory.results, name.dataset, "_clusters_on_HR_grid", "_nComp", nComp_opt, ".RData", sep = ""))

# Processed data
cluster_labels_HR <- cluster_labels_HR
cluster_labels <- cluster_labels
ARI <- ARI

# Stats
if(VERBOSE){
  cat("\nStats")
  cat(paste("\n- ARI: ", ARI))
}

# Plot cluster HR
if(PLOT){
  plot <- plot.field_points(grid, cluster_labels_HR, colormap = "H", size = 1, discrete = TRUE) +
    ggtitle("HR Clustering") + xlab("") + ylab("")
  ggsave(paste(directory.images, "clusters_HR", "_nComp", nComp_opt, ".jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}

# Plot cluster
if(PLOT){
  plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 4, discrete = TRUE) +
    ggtitle("Clustering") + xlab("") + ylab("")
  ggsave(paste(directory.images, "clusters", "_nComp", nComp_opt, ".jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Analysis DLPFC dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
PLOT <- TRUE
fdaPDE <- TRUE

setwd("~/ST-lens")
source("src/cat_utilities.R")

cat.script_title("Analysis DLPFC dataset")


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
directory.initial_data <- "data/DLPFC/"
directory.results <- "results/DLPFC/"
directory.images <- "images/DLPFC/"
name.dataset <- "DLPFC_sample9"

# Code flow control
RUN <- list()
RUN[["Pre-Processing"]] <- FALSE
RUN[["Mesh Generation"]] <- FALSE
RUN[["Mean estimation"]] <- FALSE
RUN[["Optimal nComp selection"]] <- FALSE
RUN[["Clustering at locations"]] <- FALSE


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

# Removing NaN
any(is.na(true_labels$true_label))
names.locations <- names.locations.initial[!is.na(true_labels$true_label)]
locations <- locations[names.locations,]
counts <- counts[, names.locations]
true_labels <- data.frame(true_label = true_labels[names.locations, ])
rownames(true_labels) <- names.locations

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
h <- 6.4
bbox <- NULL
seed_point <- SpatialPoints(data.frame(x = 235, y = -126))
type = "hexagonal"
simplification <- 0.25
remove_holes <- FALSE
minimum_area_hole <- NULL
simplification_hole <- 0.3
maximum_area <- 120

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
         plot = plot, width = 5, height = 6, dpi = 200)
}

# Plot mesh
if(PLOT){
  plot <- plot.fdaPDE_mesh(mesh)
  plot <- plot + xlab("") + ylab("") + ggtitle("Mesh")
  ggsave(paste(directory.images, "mesh.jpg", sep = ""),
         plot = plot, width = 5, height = 6, dpi = 200)
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
true_labels <- data.frame(true_label = true_labels[names.locations, ])
row.names(true_labels) <- names.locations

# HR grid
grid <- square_grid(SpatialPoints(locations)@bbox, 3, seed_point = seed_point)
grid <- grid[!is.na(over(SpatialPoints(grid), lattice$domain)),]

# Plot
if(PLOT){
  plot <- ggplot() +
    standard_plot_settings() + 
    xlab("x") + ylab("y") + ggtitle("HR grid") +
    geom_sf(data = st_as_sf(SpatialPoints(grid)), color = "black", size = 0.25) +
    geom_sf(data = st_as_sf(SpatialPoints(locations)), color = "red", size = 0.25)
  ggsave(paste(directory.images, "HR_grid.jpg", sep = ""),
         plot = plot, width = 5, height = 6, dpi = 200)
}


## Mean estimation ----
## ||||||||||||||||||||

cat.subsection_title("Mean estimation")

# Hyper-parameters
lambdas <- 10^seq(-4, 3, by = 0.2)

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
         plot = plot, width = 5, height = 6, dpi = 200)
}


## Optimal nComp selection ----
## ||||||||||||||||||||||||||||

# Hyper-parameters
lambda <- 10^1.8 # lambda_opt
nComp <- 20
nComp_opt <- 8

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

# Plot components HR selected
if(PLOT){
  loadings_HR <- NULL
  for(h in 1:nComp_opt){
    loadings_HR <- cbind(loadings_HR, field.eval(grid, loadings_nodes[,h], mesh))
  }
  plot <- plot.components(grid, loadings_HR, type = "tile", ncol = nComp_opt)
  ggsave(paste(directory.images, "components_selected_HR.jpg", sep = ""),
         plot = plot, width = 3*nComp_opt, height = 3*1, dpi = 200)
}

# Clean
rm(residuals_norm, loadings_nodes)


# Clustering ----
# |||||||||||||||

# Hyper-parameters
nComp_opt <- 16 # nComp_opt
clusternum <- 7
knearest <- 220 # round(sqrt(nrow(locations)))

# Clustering at locations
tic()
if(RUN[["Clustering at locations"]]){
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
  # Assembling data for clustering
  data.clustering <- NULL
  for(i in 1:nComp){
    data.clustering <- cbind(data.clustering, loadings[,i])
  }
  
  # Clustering
  cluster_labels <- walktrap_clustering(clusternum = clusternum,
                                        latent_dat = t(data.clustering[,1:nComp_opt]),
                                        knearest = knearest)
  cluster_labels_refined <- refine_cluster_10x(cluster_labels,
                                               locations,
                                               shape = "hexagon")
  
  
  # Performance index
  ARI <- adjustedRandIndex(true_labels$true_label,
                           cluster_labels)
  names(cluster_labels_refined) <- names.locations
  ARI_refiend <- adjustedRandIndex(true_labels$true_label,
                                   cluster_labels_refined)
  
  # Save clusters
  save(cluster_labels, cluster_labels_refined,
       ARI, ARI_refiend,
       # Saving options
       file = paste(directory.results, name.dataset, "_clusters_at_locations", ".RData", sep = ""))
  
  
  # Clean
  rm(cluster_labels, cluster_labels_refined,
     ARI, ARI_refiend,
     names.locations_not_unknown)
  
  # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
}
elapsed <- toc(quiet = TRUE)
times[["Clustering at locations"]] <- as.numeric(elapsed$toc - elapsed$tic)

# Load generated data
load(paste(directory.results, name.dataset, "_clusters_at_locations.RData", sep = ""))

# Processed data
cluster_labels <- cluster_labels
cluster_labels_refined <- cluster_labels_refined
ARI <- ARI
ARI_refiend <- ARI_refiend

# Stats
if(VERBOSE){
  cat("\nStats")
  cat(paste("\n- ARI:", ARI))
  cat(paste("\n- ARI refined: ", ARI_refiend))
}

# Plot cluster
if(PLOT){
  plot <- plot.field_points(locations, cluster_labels, colormap = "H", size = 1, discrete = TRUE) +
    ggtitle("Clustering") + xlab("") + ylab("")
  ggsave(paste(directory.images, "clusters_locations.jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}

# Plot cluster refined
if(PLOT){
  plot <- plot.field_points(locations, cluster_labels_refined, colormap = "H", size = 1, discrete = TRUE) +
    ggtitle("Clustering refined") + xlab("") + ylab("")
  ggsave(paste(directory.images, "clusters_refined_locations.jpg", sep = ""),
         plot = plot, width = 5, height = 5, dpi = 200)
}

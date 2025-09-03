# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Analysis HER2 dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

## ||||||||||||||||||||
# libraries ----
## ||||||||||||||||||||

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

# ARI performance index
library(mclust)

# Strings
library(glue)




#=============================================================================#




## ||||||||||||||||||||
# functions ----
## ||||||||||||||||||||

source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/plots.R")




#=============================================================================#




## ||||||||||||||||||||
# global variables ----
## ||||||||||||||||||||

cat.section_title("Global variables")

## dataset name
name_dataset <- "HER2"

## directories
path_data <- paste("data/", name_dataset, "/", sep = "")
path_images <- paste("images/", name_dataset, "/", sep = "")
mkdir(c(path_images))
path_images <- paste("images/", name_dataset, "/analysis/", sep = "")
mkdir(c(path_images))
path_results <- paste("results/", name_dataset, "/", sep = "")
mkdir(c(path_results))


## figures dimensions
figure_width <- 8
figure_height <- 8



#=============================================================================#




## ||||||||||||||||||||
# data ----
## ||||||||||||||||||||

## counts and locations ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## load
cat.section_title("Data")
load(paste(path_data, "preprocessed_data.RData", sep = ""))

# data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations
# - true_labels:    (data.frame)    n_locations x 1  (NULL if not available)


## mesh ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## load
cat.section_title("Mesh")
load(paste(path_data, "mesh.RData", sep = ""))

## stats
cat("\nStats")
cat(paste("\n- Number of locations:", nrow(locations)))
cat(paste("\n- Number of elements: ", nrow(mesh$triangles)))
cat(paste("\n- Number of nodes: ", nrow(mesh$nodes)))

## plot mesh
plot <- plot.fdaPDE_mesh(mesh) + ggtitle("Mesh")
plot


## data/mesh alignment ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## removing locations outside the domain
indexes.discarded_locations <- is.na(over(SpatialPoints(locations), lattice$domain))
names_locations <- names_locations[!indexes.discarded_locations]
locations <- locations_initial[names_locations,]
counts <- counts[, names_locations]
true_labels <- data.frame(true_label = true_labels_initial[names_locations, ])
row.names(true_labels) <- names_locations

## stats
cat("\nStats")
cat(paste("\n- Initial number of locations:", nrow(locations_initial)))
cat(paste("\n- Final number of locations:", nrow(locations)))
cat(paste("\n- Initial number of genes: ", nrow(counts_initial)))
cat(paste("\n- Final number of genes: ", nrow(counts)))

## plot final locations
plot <- plot.final_locations(SpatialPoints(locations_initial),
                             SpatialPoints(locations), lattice) +
  ggtitle("Final locations")
plot




#=============================================================================#




## ||||||||||||||||||||
# fPCA analysis ----
## ||||||||||||||||||||


## Prepare data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat.section_title("fPCA analysis")

## Save input data for time trials script
filepath <- paste(path_results, "HER2__domain_input_for_FPCA.RData", sep = "")
save(mesh, grid, counts, locations, file=filepath)

## HR grid
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
       plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")

## prepare data
mesh_data <- list(
  nodes = as.matrix(mesh$nodes),
  edges = as.matrix(mesh$edges),
  elements = as.matrix(mesh$triangles),
  neigh = as.matrix(mesh$neighbors),
  boundary = as.matrix(as.numeric(mesh$nodesmarkers))
)
domain <- fdaPDE2::Mesh(mesh_data)
data <- fdaPDE2::functional_data(
  domain = domain,
  X = counts, 
  locations = locations
)





## Run fPCA, GCV ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## start timer
start.time <- Sys.time()

## model parameters
n_comp <- 20

## lambda selected
lambda_grid <-  hyperparameters(10^seq(-9, 2, by = 0.2))

## fPCA model initialization
model_fPCA_gcv <- fdaPDE2::fPCA(
  data,
  center = centering(calibrator = gcv(lambda_grid, seed=0)),
  solver = sequential()
)

## fPCA model fit
model_fPCA_gcv$fit(
  calibrator = gcv(lambda_grid, seed = 0),
  n_pc = n_comp
)

## stop timer
end.time <- Sys.time()
execution_time_gcv <- end.time - start.time
cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))





## Plot fPCA: GCV ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Set model 
model <- model_fPCA_gcv


## mean
mean <- model$results$X_mean
mean_locs <- evaluate_field(locations, mean, mesh)
mean_HR <- evaluate_field(grid, mean, mesh)

## plot mean tile
plot <- plot.field_tile_original(grid, mean_HR, LEGEND = TRUE) +
  standard_plot_settings_fields() +
  ggtitle("Mean: interpolated")
plot
ggsave(paste(path_images, "2_mean_GCV.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## plot mean at locations
plot <- plot.field_points(locations, mean_locs, size = 2) +
  standard_plot_settings_fields() +
  ggtitle("Mean: at locations")
plot
ggsave(paste(path_images, "2_mean_at_locations_GCV.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## loadings re-normalization
loadings_NotNormalized <- model$results$loadings
scores_NotNormalized <- model$results$scores
loadings <- model$results$loadings
scores <- model$results$scores
loadings_locs <- NULL
loadings_HR <- NULL
loadings_locs_NotNormalized <- NULL
loadings_HR_NotNormalized <- NULL
changes_of_sign <- rep(1, n_comp)
for(i in 1:n_comp) {
  loadings_HR_NotNormalized <- cbind(loadings_HR_NotNormalized, evaluate_field(grid, loadings_NotNormalized[, i], mesh))
  loadings_locs_NotNormalized <- cbind(loadings_locs_NotNormalized, evaluate_field(locations, loadings_NotNormalized[, i], mesh))
  loadings[, i] <- changes_of_sign[i] * loadings[, i] * norm_l2(scores[, i])
  scores[, i] <- scores[, i]/norm_l2(scores[, i])
  loadings_HR <- cbind(loadings_HR, evaluate_field(grid, loadings[, i], mesh))
  loadings_locs <- cbind(loadings_locs, evaluate_field(locations, loadings[ ,i], mesh))
}


## plot components tile (scaled by limits of ALL loadings)
n_comp_plot <- 4
plot_list <- list()
limits = range(loadings_HR)
for(i in 1:n_comp_plot) {
  plot_list[[i]] <- plot.field_tile(grid, loadings_HR[, i], limits = limits) +
    standard_plot_settings_fields() +
    ggtitle(paste("w", i, sep = "")) +
    #theme(legend.position = "none") +
    NULL
}
plot <- arrangeGrob(grobs = plot_list, nrow = ceiling(sqrt(n_comp_plot)))
grid.arrange(plot)
ggsave(paste(path_images, glue("{n_comp_plot}_components_GCV.pdf"), sep = ""),
       plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")

## plot components at locations (scaled by limits of ALL loadings)
plot_list <- list()
limits = range(loadings_locs)
for(i in 1:n_comp_plot) {
  plot_list[[i]] <- plot.field_points(locations, loadings_locs[, i], limits = limits, size = 1.7) +
    standard_plot_settings_fields() +
    ggtitle(paste("w", i, sep = "")) +
    #theme(legend.position = "none") +
    NULL
}
plot <- arrangeGrob(grobs = plot_list, nrow = ceiling(sqrt(n_comp_plot)))
grid.arrange(plot)
ggsave(paste(path_images, glue("{n_comp_plot}_components_at_locations_GCV.pdf"), sep = ""),
       plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")



## Save fPCA +GCV results ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_gcv <- mean
scores_gcv <- scores
loadings_gcv <- loadings
loadings_HR_gcv <- loadings_HR
loadings_locs_gcv <- loadings_locs
loadings_locs_gcv_NotNormalized <- loadings_locs_NotNormalized
loadings_HR_gcv_NotNormalized <- loadings_HR_NotNormalized
save(
  execution_time_gcv, n_comp,
  mean_gcv, scores_gcv, 
  loadings_gcv, loadings_HR_gcv, loadings_locs_gcv,
  loadings_locs_gcv_NotNormalized, loadings_HR_NotNormalized,
  file = paste(path_results, "fPCA_gcv.RData", sep = "")
)





#=============================================================================#



## ||||||||||||||||||||
# saving non-fPCA results ----
## ||||||||||||||||||||

## final data
save(
  ## initial data
  locations_initial,   names_locations_initial,
  counts_initial, names_genes.initial,
  true_labels_initial,
  ## analyzed data
  locations, names_locations,
  counts, names_genes,
  true_labels,
  ## file name
  file = paste(path_data, "/analyzed_data.RData", sep = "")
)



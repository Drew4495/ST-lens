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

# ARI performance index
library(mclust)

# functions ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/plots.R")


##############################################################################


cat.script_title("Analysis HER2 dataset")


# global variables ----

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


# data ----

## counts and locations ----

## load
cat.section_title("Data")
load(paste(path_data, "preprocessed_data.RData", sep = ""))

# data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations
# - true_labels:    (data.frame)    n_locations x 1  (NULL if not available)


## mesh ----

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


# fPCA analysis ----

cat.section_title("fPCA analysis")

## HR grid
grid_step <- 1/6
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



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## start timer
start.time <- Sys.time()

## model parameters
n_comp <- 4

## lambda selected
lambda_selected <- hyperparameters(1e-9)

## fPCA model initialization
model_fPCA <- fdaPDE2::fPCA(
  data,
  center = centering(lambda_selected),
  solver = sequential()
)

## fPCA model fit
model_fPCA$fit(
  calibrator = lambda_selected,
  n_pc = n_comp
)

## stop timer
end.time <- Sys.time()
execution_time <- end.time - start.time
cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## mean
mean <- model_fPCA$results$X_mean
mean_locs <- evaluate_field(locations, mean, mesh)
mean_HR <- evaluate_field(grid, mean, mesh)

## plot mean tile
plot <- plot.field_tile(grid, mean_HR, LEGEND = TRUE) +
  standard_plot_settings_fields() +
  ggtitle("Mean")
plot
ggsave(paste(path_images, "2_mean.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## plot mean at locations
plot <- plot.field_points(locations, mean_locs, size = 2) +
  standard_plot_settings_fields() +
  ggtitle("Mean")
plot
ggsave(paste(path_images, "2_mean_at_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## loadings re-normalization
loadings <- model_fPCA$results$loadings
scores <- model_fPCA$results$scores
loadings_locs <- NULL
loadings_HR <- NULL
changes_of_sign <- c(1,1,1,1)
for(i in 1:n_comp) {
  loadings[, i] <- changes_of_sign[i] * loadings[, i] * norm_l2(scores[, i])
  scores[, i] <- scores[, i]/norm_l2(scores[, i])
  loadings_HR <- cbind(loadings_HR, evaluate_field(grid, loadings[, i], mesh))
  loadings_locs <- cbind(loadings_locs, evaluate_field(locations, loadings[, i], mesh))
}

## plot components tile
plot_list <- list()
limits = range(loadings_HR)
for(i in 1:n_comp) {
  plot_list[[i]] <- plot.field_tile(grid, loadings_HR[, i], limits = limits) +
    standard_plot_settings_fields() +
    ggtitle(paste("w", i, sep = ""))
}
plot <- arrangeGrob(grobs = plot_list, nrow = ceiling(sqrt(n_comp)))
grid.arrange(plot)
ggsave(paste(path_images, "3_components.pdf", sep = ""),
       plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")


## plot components at locations
plot_list <- list()
limits = range(loadings_locs)
for(i in 1:n_comp) {
  plot_list[[i]] <- plot.field_points(locations, loadings_locs[, i], limits = limits, size = 3) +
    standard_plot_settings_fields() +
    ggtitle(paste("w", i, sep = ""))
}
plot <- arrangeGrob(grobs = plot_list, nrow = ceiling(sqrt(n_comp)))
grid.arrange(plot)
ggsave(paste(path_images, "3_components_at_locations.pdf", sep = ""),
       plot = plot, width = figure_width*2, height = figure_height*2, dpi = "print", units = "cm")


# saving results ----

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

## analysis results
save(
  execution_time, n_comp,
  mean, scores, loadings,
  file = paste(path_results, "fPCA.RData", sep = "")
)

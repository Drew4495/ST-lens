# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Example: HER2 mesh generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
fdaPDE <- FALSE

# setwd("~/ST-lens")
source("src/cat_utilities.R")

cat.script_title("Example: HER2 mesh generation")


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


# |||||||||||||||||||||
# Global variables ----
# |||||||||||||||||||||

cat.section_title("Global variables")

# Images directory
directory.images <- "images/example_mesh_generation/"
if (!file.exists(directory.images)){
  dir.create(directory.images)
}
directory.images <- paste(directory.images, "/HER2/", sep = "")
if (!file.exists(directory.images)){
  dir.create(directory.images)
}


# |||||||||
# Data ----
# |||||||||

load("data/HER2/HER2.RData")

locations <- SpatialPoints(locations)

# range(locations@coords[,1])
xmin <- 0
xmax <- 35 
# range(locations@coords[,2])
ymin <- -36
ymax <- -7

# Plot
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(directory.images, "0_initial_locations.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")


# ||||||||||||||||||||||||||
# Domain identification ----
# ||||||||||||||||||||||||||

## Grid ----
## |||||||||

# Lattice type
type = "square"

# Step of the grid
h <- 1
# It is decided by the user by looking at the initial distribution of location. 
# Lower the value of h more precise will be the domain reconstruction.
# However, it can not be too low otherwise there could be unwelcome holes in 
# the domain

# Seed Point
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))
# It is the seed for the generation of the grid, the final grid is guaranteed to
# contain this point. It is used to have always the same grid.

# Plot
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 0.15) +
  geom_sf(data = st_as_sf(square(seed_point, h/sqrt(2))), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "1_initial_locations_check.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")

# Grid generation
grid <- generate_grid(locations@bbox, h, seed_point, type = type)

# Plot
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(grid), color = "blue", size = 0.05) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 0.1)
ggsave(paste(directory.images, "2_grid.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")


## Lattice ----
## ||||||||||||

check <- c()
polygons <- list()
polygons_all <- list()
for(i in 1:nrow(grid@coords)){
  point <- grid[i,]
  polygon <- square(point, h/sqrt(2) - 1e-9)
  check[i] <- any(!is.na(over(locations, polygon)))
  if(check[i]){
    polygons <- c(polygons, polygon@polygons[[1]]@Polygons[[1]])
  }
  polygons_all <- c(polygons_all, polygon@polygons[[1]]@Polygons[[1]])
}

lattice <- SpatialPolygons(list(Polygons(polygons, ID = "lattice")))
lattice_all <- SpatialPolygons(list(Polygons(polygons_all, ID = "lattice")))

# Plot all
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "3_lattice.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")

# Plot all and locations
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black") +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(directory.images, "4_lattice_and_locations.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")

# Plot selected
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice), fill = "blue", alpha = 0.5, color = "black")+
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(directory.images, "5_lattice_only_selected.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")


## Simplification ----
## |||||||||||||||||||

# User defiend parameter
simplification <- 0.15
# This parameter represents the percentage of boundary vertices to be kept.
# The user should set a value such that the boundary is simplified enough
# but without exeeding otherwise there will be a lot of discarded points

# About holes
remove_holes <- FALSE
minimum_area_hole <- NULL
simplification_hole <- 1

# Lattice
lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = "square")

# Plot
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "6_domain.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")

# Simplification
lattice_simplified <- simplify_domain(lattice, simplification,
                                      remove_holes, minimum_area_hole, simplification_hole)

# During this step the islands, namely the parts of domain that are not joined 
# to the one with largest, area are discarded

# Plot
plot <- plot.discretized_domain(lattice_simplified, size = 0.1) +
  mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + guides(color = "none")
ggsave(paste(directory.images, "7_domain_simplified.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")

# Discarded points
indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
locations.final <- locations[!indexes.discarded_locations,]

# Plot
plot <- plot.final_locations(locations, SpatialPoints(locations.final), lattice_simplified, size = 0.05) +
  mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(directory.images, "8_final_locations.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")


# ||||||||||||
# Meshing ----
# ||||||||||||

# User defiend parameter
maximum_area <- 0.3
# It is the threshold for the larger possible value for an element of the mesh.
# Is the original mesh contain elements larger than it it is refined until all 
# the elements meet this constraint.
# The user should set a value such that the final number of nodes of the mesh
# has the same order of magnitude of the number of locations.

# Mesh generation
mesh <- generate_mesh(lattice_simplified, maximum_area)

# Plot
plot <- plot.fdaPDE_mesh(mesh) + 
  mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) 
ggsave(paste(directory.images, "9_mesh.pdf", sep = ""),
       plot = plot, width = 6, height = 5, dpi = "print", units = "cm")



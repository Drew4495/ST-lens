# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Example: DLPFC mesh generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
fdaPDE <- FALSE

# setwd("~/ST-lens")
source("src/cat_utilities.R")

cat.script_title("Example: DLPFC mesh generation")


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
directory.images <- paste(directory.images, "/DLPFC_noNaN/", sep = "")
if (!file.exists(directory.images)){
  dir.create(directory.images)
}


# |||||||||
# Data ----
# |||||||||

load("data/DLPFC/DLPFC_sample9.RData")

# Removing NaN
any(is.na(true_labels$true_label))
locations <- locations[!is.na(true_labels$true_label),]
locations <- SpatialPoints(locations)

# range(locations@coords[,1])
xmin <- 125
xmax <- 510 
# range(locations@coords[,2])
ymin <- -515
ymax <- -110 

# Plot
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Initial locations") + 
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 1)
ggsave(paste(directory.images, "0_initial_locations.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)


# ||||||||||||||||||||||||||
# Domain identification ----
# ||||||||||||||||||||||||||

## Grid ----
## |||||||||

# Lattice type
type = "hexagonal"

# Why an hexagonal grid?
# https://strimas.com/post/hexagonal-grids/

# Regular hexagons are the closest shape to a circle that can be used for the
# regular tessellation of a plane and they have additional symmetries compared 
# to squares. These properties give rise to the following benefits.

# - Reduced edge effects: a hexagonal grid gives the lowest perimeter to area 
#   ratio of any regular tessellation of the plane. In practice, this means that 
#   edge effects are minimized when working with hexagonal grids. This is 
#   essentially the same reason beehives are built from hexagonal honeycomb: it 
#   is the arrangement that minimizes the amount of material used to create a 
#   lattice of cells with a given volume.
# - All neighbours are identical: square grids have two classes of neighbours, 
#   those in the cardinal directions that share an edge and those in diagonal 
#   directions that share a vertex. In contrast, a hexagonal grid cell has six 
#   identical neighbouring cells, each sharing one of the six equal length 
#   sides. Furthermore, the distance between centroids is the same for all 
#   neighbours.
# - Better fit to curved surfaces: when dealing with large areas, where the 
#   curvature of the earth becomes important, hexagons are better able to fit 
#   this curvature than squares. This is why soccer balls are constructed of 
#   hexagonal panels.
# - They look badass: it can’t be denied that hexagonal grids look way more 
#   impressive than square grids!

# Step of the grid
h <- 6.4
# It is decided by the user by looking at the initial distribution of location. 
# Lower the value of h more precise will be the domain reconstruction.
# However, it can not be too low otherwise there could be unwelcome holes in 
# the domain

# Seed Point
seed_point <- SpatialPoints(data.frame(x = 235, y = -126))
# It is the seed for the generation of the grid, the final grid is guaranteed to
# contain this point. It is used to have always the same grid.

# Plot
plot <- ggplot() +
  standard_plot_settings() + 
  xlab("x") + ylab("y") + ggtitle("Check Parameters") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 1) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 1.5) +
  geom_sf(data = st_as_sf(hex(seed_point, h/sqrt(3))), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "1_initial_locations_check.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)

# Grid generation
grid <- generate_grid(locations@bbox, h, seed_point, type = type)

# Plot
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Hexagonal Grid") + 
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(grid), color = "blue", size = 0.5) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 1)
ggsave(paste(directory.images, "2_grid.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)


## Lattice ----
## ||||||||||||

check <- c()
polygons <- list()
polygons_all <- list()
for(i in 1:nrow(grid@coords)){
  point <- grid[i,]
  polygon <- hex(point, h/sqrt(3) - 1e-9)
  check[i] <- any(!is.na(over(locations, polygon)))
  if(check[i]){
    polygons <- c(polygons, polygon@polygons[[1]]@Polygons[[1]])
  }
  polygons_all <- c(polygons_all, polygon@polygons[[1]]@Polygons[[1]])
}

lattice <- SpatialPolygons(list(Polygons(polygons, ID = "hex_lattice")))
lattice_all <- SpatialPolygons(list(Polygons(polygons_all, ID = "hex_lattice")))

# Plot all
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Lattice") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "3_lattice.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)

# Plot all and locations
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Lattice and locations") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black") +
  geom_sf(data = st_as_sf(locations), color = "black", size = 1)
ggsave(paste(directory.images, "4_lattice_and_locations.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)

# Plot selected
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Selected hexagons") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "5_lattice_only_selected.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)


## Simplification ----
## |||||||||||||||||||

# User defiend parameter
simplification <- 0.25
# This parameter represents the percentage of boundary vertices to be kept.
# The user should set a value such that the boundary is simplified enough
# but without exeeding otherwise there will be a lot of discarded points

# About holes
remove_holes <- FALSE
minimum_area_hole <- NULL
simplification_hole <- 0.3

# Lattice
lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = type)

# Plot
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Domain") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "6_domain.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)

# Simplification
lattice_simplified <- simplify_domain(lattice, simplification,
                                      remove_holes, minimum_area_hole, simplification_hole)

# During this step the islands, namely the parts of domain that are not joined 
# to the one with largest, area are discarded

# Plot
plot <- plot.discretized_domain(lattice_simplified) +
  xlab("") + ylab("") + ggtitle("Domain simplified") +
  xlim(xmin, xmax) + ylim(ymin, ymax) + guides(color = "none")
ggsave(paste(directory.images, "7_domain_simplified.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)

# Discarded points
indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
locations.final <- locations[!indexes.discarded_locations,]

# Plot
plot <- plot.final_locations(locations, SpatialPoints(locations.final), lattice_simplified, size = 0.5) +
  xlab("") + ylab("") + ggtitle("Final locations") +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(directory.images, "8_final_locations.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)


# ||||||||||||
# Meshing ----
# ||||||||||||

# User defiend parameter
maximum_area <- 120
# It is the threshold for the larger possible value for an element of the mesh.
# Is the original mesh contain elements larger than it it is refined until all 
# the elements meet this constraint.
# The user should set a value such that the final number of nodes of the mesh
# has the same order of magnitude of the number of locations.

# Mesh generation
mesh <- generate_mesh(lattice_simplified, maximum_area)

# Plot
plot <- plot.fdaPDE_mesh(mesh) + 
  xlab("") + ylab("") + ggtitle("Mesh") +
  xlim(xmin, xmax) + ylim(ymin, ymax) 
ggsave(paste(directory.images, "9_mesh.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)



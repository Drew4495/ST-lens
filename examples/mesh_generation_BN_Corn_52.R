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
directory.images <- paste(directory.images, "/BN_Corn_52/", sep = "")
if (!file.exists(directory.images)){
  dir.create(directory.images)
}


# |||||||||
# Data ----
# |||||||||

filepath.data.genexcell <- "data/raw/cell_segmentation/BN-corn-52_HighThreshold_gene_x_cell.csv"
filepath.data.cell_locs <- "data/raw/cell_segmentation/BN-corn-52_HighThreshold_cells.csv"
#Load in counts
counts.initial <- read.csv(filepath.data.genexcell) 
counts_colnames <- rownames(counts.initial[-1,-1])
counts.initial <- subset(counts.initial, select= -c(cell_x, cell_y, cell_area))
counts.initial <- t(apply(counts.initial[-1, -1], 2, as.integer)) #cell_id = 0 are spots not assigned to cells
colnames(counts.initial) <- counts_colnames

#Load in locations
locations.cells <- read.csv(filepath.data.cell_locs)
locations.initial <- subset(locations.cells, select=-area)
rm(locations.cells)
rownames(locations.initial) <- locations.initial$cell_id
locations.initial <- subset(locations.initial, select=-cell_id)


#locations <- SpatialPoints(locations.initial)
locations <- SpatialPoints(locations.significant)

# range(locations@coords[,1])
xmin <- 800
xmax <- 16786
# range(locations@coords[,2])
ymin <- 1458
ymax <- 32309

# Plot
plot <- ggplot() + mesh_generation_examples_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(directory.images, "0_initial_locations.pdf", sep = ""),
       plot = plot, width = 7.5, height = 8, dpi = "print", units = "cm")


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
h <- 250
# It is decided by the user by looking at the initial distribution of location. 
# Lower the value of h more precise will be the domain reconstruction.
# However, it can not be too low otherwise there could be unwelcome holes in 
# the domain

# Seed Point
seed_point <- SpatialPoints(data.frame(x = 10000, y = 20000))
# It is the seed for the generation of the grid, the final grid is guaranteed to
# contain this point. It is used to have always the same grid.

# Plot
plot <- ggplot() + #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 0.15) +
  geom_sf(data = st_as_sf(hex(seed_point, h/sqrt(3))), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "1_initial_locations_check.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")

# Grid generation
grid <- generate_grid(locations@bbox, h, seed_point, type = type)

# Plot
plot <- ggplot() + #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(grid), color = "blue", size = 0.05) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 1)
ggsave(paste(directory.images, "2_grid.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")


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
plot <- ggplot() + #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "3_lattice.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")

# Plot all and locations
plot <- ggplot() + #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black") +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(directory.images, "4_lattice_and_locations.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")

# Plot selected
plot <- ggplot() + #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "5_lattice_only_selected.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")


## Simplification ----
## |||||||||||||||||||

# User defiend parameter
simplification <- 0.25
# This parameter represents the percentage of boundary vertices to be kept.
# The user should set a value such that the boundary is simplified enough
# but without exeeding otherwise there will be a lot of discarded points

# About holes
remove_holes <- FALSE
minimum_area_hole <- 200000
simplification_hole <- 0.5

# Lattice
lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = type)

# Plot
plot <- ggplot() + #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "6_domain.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")

# Simplification
lattice_simplified <- simplify_domain(lattice, simplification,
                                      remove_holes, minimum_area_hole, simplification_hole)

# During this step the islands, namely the parts of domain that are not joined 
# to the one with largest, area are discarded

# Plot
plot <- plot.discretized_domain(lattice_simplified, size = 0.1) +
  #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + guides(color = "none")
ggsave(paste(directory.images, "7_domain_simplified.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")

# Discarded points
indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
locations.final <- locations[!indexes.discarded_locations,]

# Plot
plot <- plot.final_locations(locations, SpatialPoints(locations.final), lattice_simplified, size = 0.05) +
  #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(directory.images, "8_final_locations.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")


# ||||||||||||
# Meshing ----
# ||||||||||||

# User defiend parameter
maximum_area <- 100000
# It is the threshold for the larger possible value for an element of the mesh.
# Is the original mesh contain elements larger than it it is refined until all 
# the elements meet this constraint.
# The user should set a value such that the final number of nodes of the mesh
# has the same order of magnitude of the number of locations.

# Mesh generation
mesh <- generate_mesh(lattice_simplified, maximum_area)

# Plot
plot <- plot.fdaPDE_mesh(mesh) + 
  #mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) 
ggsave(paste(directory.images, "9_mesh.pdf", sep = ""),
       plot = plot,width = 7.5, height = 8, dpi = "print", units = "cm")



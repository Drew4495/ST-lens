# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Example: HER2 mesh generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE

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
plot <- ggplot() +
  standard_plot_settings() + 
  xlab("x") + ylab("y") + ggtitle("Check Parameters") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 1) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 1.5) +
  geom_sf(data = st_as_sf(square(seed_point, h/sqrt(2))), fill = "blue", alpha = 0.5, color = "black")
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
  geom_sf(data = st_as_sf(lattice), fill = "blue", alpha = 0.5, color = "black")+
  geom_sf(data = st_as_sf(locations), color = "black", size = 1)
ggsave(paste(directory.images, "5_lattice_only_selected.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)


## Simplification ----
## |||||||||||||||||||

# User defiend parameter
simplification <- 0.15
# This parameter represents the percentage of boundary vertices to be kept.
# The user should set a value such that the boundary is simplified enough
# but without exeeding otherwise there will be a lot of discarded points

# Lattice
lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = "square")

# Plot
plot <- ggplot() + standard_plot_settings() +
  xlab("") + ylab("") + ggtitle("Domain") +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(directory.images, "6_domain.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)

# Simplification
lattice_simplified <- simplify_domain(lattice, simplification)

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
  xlab("") + ylab("") + ggtitle("Mesh") +
  xlim(xmin, xmax) + ylim(ymin, ymax) 
ggsave(paste(directory.images, "9_mesh.jpg", sep = ""),
       plot = plot, width = 6.5, height = 7, dpi = 200)



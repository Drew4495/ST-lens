# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% HER2 dataset mesh generation - old version %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

# libraries ----
suppressMessages(library(sf))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(raster))
suppressMessages(library(rmapshaper))
suppressMessages(library(fdaPDE))

# sources ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/plots.R")
source("src/meshing.R")
source("src/plots.R")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cat.script_title("HER2 dataset mesh generation - old version")


# global variables ----

cat.section_title("Global variables")

name_dataset <- "HER2"

## paths
path_data <- paste("data/", name_dataset, "/", sep = "")
path_images <- paste("images/", name_dataset, "_old/", sep = "")
mkdir(c(path_images))
path_images <- paste("images/", name_dataset, "_old/meshing/", sep = "")
mkdir(c(path_images))

## figures dimensions
figure_width <- 8
figure_height <- 8


# data ----

cat.section_title("Data")

## load data
load(paste(path_data, "data.RData", sep = ""))
rm(counts, true_labels)

## create sp object
locations <- SpatialPoints(locations)

# range(locations@coords[,1])
xmin <- 1
xmax <- 34
# range(locations@coords[,2])
ymin <- -36
ymax <- -5

## plot
plot <- ggplot() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + mesh_generation_plot_settings() + 
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(path_images, "0_initial_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


# domain identification ----

cat.section_title("Domain identification")


## grid ----

cat.subsection_title("Grid generation")

## lattice type
type = "square"

## step of the grid
h <- 1

## grid seed Point
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))

## plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 0.15) +
  geom_sf(data = st_as_sf(square(seed_point, h/sqrt(3))), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "1_initial_locations_check.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## grid generation
grid <- generate_grid(locations@bbox, h, seed_point, type = type)

## plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(grid), color = "blue", size = 0.05) +
  geom_sf(data = st_as_sf(seed_point), color = "red", size = 1)
ggsave(paste(path_images, "2_grid.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## lattice ----

cat.subsection_title("Lattice generation")

## select polygons to keep
check <- c()
polygons <- list()
polygons_all <- list()
for(i in 1:nrow(grid@coords)){
  point <- grid[i,]
  polygon <- square(point, h/sqrt(3) - 1e-9)
  check[i] <- any(!is.na(over(locations, polygon)))
  if(check[i]){
    polygons <- c(polygons, polygon@polygons[[1]]@Polygons[[1]])
  }
  polygons_all <- c(polygons_all, polygon@polygons[[1]]@Polygons[[1]])
}
lattice <- SpatialPolygons(list(Polygons(polygons, ID = "hex_lattice")))
lattice_all <- SpatialPolygons(list(Polygons(polygons_all, ID = "hex_lattice")))

## plot all
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "3_lattice.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## plot all and locations
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice_all), fill = "blue", alpha = 0.5, color = "black") +
  geom_sf(data = st_as_sf(locations), color = "black", size = 0.1)
ggsave(paste(path_images, "4_lattice_and_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## plot selected
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice), fill = "blue", alpha = 0.5, color = "black")
ggsave(paste(path_images, "5_lattice_only_selected.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


## simplification ----
## |||||||||||||||||||

cat.subsection_title("Domain simplification")

## percentage of boundary vertices to be kept
simplification <- 0.15

## holes parameters
remove_holes <- FALSE
minimum_area_hole <- NULL
simplification_hole <- 1

## lattice
lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = type)

## plot
plot <- ggplot() + mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) +
  geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black")
ggsave(paste(path_images, "6_domain.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## domain simplification
lattice_simplified <- simplify_domain(lattice, simplification,
                                      remove_holes, minimum_area_hole, simplification_hole)

## plot
plot <- plot.discretized_domain(lattice_simplified, size = 0.1) +
  mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax) + guides(color = "none")
ggsave(paste(path_images, "7_domain_simplified.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")

## discarded points
indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
locations.final <- locations[!indexes.discarded_locations,]

## plot
plot <- plot.final_locations(locations, SpatialPoints(locations.final), lattice_simplified, size = 0.05) +
  mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(path_images, "8_final_locations.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


# meshing ----

cat.section_title("Meshing")

## maximum area for a mesh element
maximum_area <- 0.3

## mesh generation
mesh <- generate_mesh(lattice_simplified, maximum_area)

## plot
plot <- plot.fdaPDE_mesh(mesh) +
  mesh_generation_plot_settings() +
  xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste(path_images, "9_mesh.pdf", sep = ""),
       plot = plot, width = figure_width, height = figure_height, dpi = "print", units = "cm")


# saving results ----
lattice <- lattice_simplified
save(mesh, lattice,
     file = paste(path_data, "mesh_old.RData", sep = ""))

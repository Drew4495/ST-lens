rm(list = ls())
graphics.off()

## ||||||||||||||||||||
# Libraries and Functions ----
## ||||||||||||||||||||


# Libraries ----
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

# Timing
library(tictoc)

## From meshing script
suppressMessages(library(sf))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(raster))
suppressMessages(library(rmapshaper))
suppressMessages(library(fdaPDE))



# Functions ----
## ||||||||||||||||||||
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/meshing.R")
source("src/plots.R")




#=============================================================================#




## ||||||||||||||||||||
# Define Globals ----
## ||||||||||||||||||||

cat.section_title("Global variables")

## dataset name
name_dataset <- "DLPFC"

## directories
path_data <- paste("data/", name_dataset, "/", sep = "")
path_images <- paste("images/", name_dataset, "/", sep = "")
mkdir(c(path_images))
path_images <- paste("images/", name_dataset, "/analysis/", sep = "")
mkdir(c(path_images))
path_results <- paste("results/", name_dataset, "/", sep = "")
mkdir(c(path_results))

mesh_input_path <- paste(path_data, "data.RData", sep = "")
fpca_input_path <- paste(path_results, "domain_input_for_FPCA.RData", sep = "")




#=============================================================================#




## ||||||||||||||||||||
# Mesh time trials: Load data ----
## ||||||||||||||||||||
load(mesh_input_path)

n_trials <- 4
trial_list <- list()

for (trial_i in 1:n_trials){
  ### Clear timing log and print iteration number
  tic.clearlog()
  cat("Iteration:", trial_i, "\n")
  
  tic("Mesh Overall time")
  
  
  
  # data
  locations <- SpatialPoints(locations)
  
  # domain identification (grid, lattice, simplification)
  type = "hexagonal" ## lattice type
  h <- 6.4 ## step of the grid
  seed_point <- SpatialPoints(data.frame(x = 235, y = -126)) ## grid seed Point
  grid <- generate_grid(locations@bbox, h, seed_point, type = type) ## grid generation
  
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
  
  ## percentage of boundary vertices to be kept
  simplification <- 0.25
  
  ## holes parameters
  remove_holes <- TRUE
  minimum_area_hole <- NULL
  simplification_hole <- NULL
  
  ## lattice
  lattice <- generate_lattice(locations, h, locations@bbox, seed_point, type = type)
  
  ## domain simplification
  lattice_simplified <- simplify_domain(lattice, simplification,
                                        remove_holes, minimum_area_hole, simplification_hole)
  
  ## discarded points
  indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
  locations.final <- locations[!indexes.discarded_locations,]
  
  # Meshing
  maximum_area <- 120 ## maximum area for a mesh element
  mesh <- generate_mesh(lattice_simplified, maximum_area) ## mesh generation
  
  toc(log=T)
  
  
  
  ### Store time trials
  times_log <- tic.log(format=F)
  trial_df <- data.frame(
    Iteration = trial_i,
    Step = sapply(times_log, function(t) t$msg),
    ElapsedTime = sapply(times_log, function(t) t$toc[[1]] - t$tic[[1]])
  )
  trial_list[[trial_i]] <- trial_df
  
}


mesh_timetrials_df <- bind_rows(trial_list)

filepath <- paste(path_results, "timetrials_mesh_df.RData", sep="")
save(mesh_timetrials_df, file=filepath)






#=============================================================================#




## ||||||||||||||||||||
# fPCA time trials: Load data ----
## ||||||||||||||||||||
load(fpca_input_path)


## Prepared data (because issues with saving/loading post-prepared data)
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



## Run trials
n_trials <- 4
trial_list <- list()

for (trial_i in 1:n_trials){
  ### Clear timing log and print iteration number
  tic.clearlog()
  cat("Iteration:", trial_i, "\n")
  
  tic("fpca: overall time")
  tic("fpca: model initialization")
  
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
  toc(log=T)
  
  tic("fpca: model fitting")
  ## fPCA model fit
  model_fPCA_gcv$fit(
    calibrator = gcv(lambda_grid, seed = 0),
    n_pc = n_comp
  )
  toc(log=T)
  
  toc(log=T)
  
  
  
  ### Store time trials
  times_log <- tic.log(format=F)
  trial_df <- data.frame(
    Iteration = trial_i,
    Step = sapply(times_log, function(t) t$msg),
    ElapsedTime = sapply(times_log, function(t) t$toc[[1]] - t$tic[[1]])
  )
  trial_list[[trial_i]] <- trial_df
  
}


fpca_timetrials_df <- bind_rows(trial_list)
filepath <- paste(path_results, "timetrials_fpca_df.RData", sep="")
save(fpca_timetrials_df, file=filepath)













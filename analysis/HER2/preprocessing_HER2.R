# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Preprocessing HER2 dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

## SpatialPCA
library(SpatialPCA)

## fdaPDE
library(fdaPDE)


# functions ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/plots.R")


##############################################################################


cat.script_title("Preprocessing HER2 dataset")


# global variables ----

cat.section_title("Global variables")

## dataset name
name_dataset <- "HER2"

## paths
path_data <- paste("data/", name_dataset, "/", sep = "")


# data ----

## load
cat.section_title("Data")
load(paste(path_data, "data.RData", sep = ""))

# Data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations
# - true_labels:    (data.frame)    n_locations x 1  (NULL if not available)

## save initial data
locations_initial <- locations
counts_initial <- counts
true_labels_initial <- true_labels
names_locations_initial <- rownames(locations)
names_genes.initial <- rownames(counts)


# preprocessing ----

cat.section_title("Preprocessing")

## preprocessing
ST = CreateSpatialPCAObject(
  counts = counts,
  location = as.matrix(locations),
  gene.type = "spatial",
  sparkversion = "spark",
  gene.number = 3000,
  customGenelist = NULL,
  min.loctions = 20,
  min.features = 20
)

## processed data
locations <- data.frame(ST@location)
counts <- ST@normalized_expr
names_locations <- rownames(locations)
names_genes <- rownames(counts)
rm(ST)

# scale the expression of each gene
for(i in 1:nrow(counts)){
  counts[i,] = scale(counts[i,])
}

# saving results ----
save(
  ## initial data
  locations_initial,   names_locations_initial,
  counts_initial, names_genes.initial,
  true_labels_initial,
  ## data after preprocessing
  locations, names_locations,
  counts, names_genes,
  true_labels,
  ## file name
  file = paste(path_data, "/preprocessed_data.RData", sep = "")
)

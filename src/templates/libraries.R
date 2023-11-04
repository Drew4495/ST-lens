# Libraries ----
# ||||||||||||||

# Statistical utilities
library(MASS)

# Spatial
library(sp)
library(sf)
library(rmapshaper)
library(spatstat)

# Plots
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)

# Time
library(tictoc)

# Pre-Processing
library(Seurat)
library(SPARK)

# Mesh
library(fdaPDE)

# Data Decomposition
if(fdaPDE){
  library(fdaPDE2)
}

# ARI performance index
library(mclust)
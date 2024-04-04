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
library(pals)
library(RColorBrewer)
library(ComplexHeatmap)

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

# ARI and CHAOS performance index
library(mclust)
library(parallel)
library(pdist)
library(lisi)

# Coding Ease
library(glue)
library(reshape2)
library(tibble)
library(dplyr)
library(tidyr)

#Stats
library(car)
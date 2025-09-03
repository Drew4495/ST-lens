# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% SlideseqCerebellum dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()





## ||||||||||||||||||||
# libraries and functions ----
## ||||||||||||||||||||

## plots
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(patchwork)

## Spatial data
library(sp)
library(sf)
library(rmapshaper)
library(spatstat)

## fdaPDE
library(fdaPDE)
library(fdaPDE2)

## ARI performance index
library(mclust)

## parallel
library(foreach)
library(doParallel)

## Other utils
library(glue)


## functions
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/geometry.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/plots.R")
source("src/clustering.R")
source("src/src_alignment.R")
source("src/src_NumComps_Selection.R")




#=============================================================================#




## ||||||||||||||||||||
# Set globals ----
## ||||||||||||||||||||

cat.script_title("Deciding Number of Spatial Components (NSC) for clustering: DLPFC dataset")

# global variables 

cat.section_title("Global variables")

## dataset name
name_dataset <- "SlideseqCerebellum"

## directories
path_data <- paste("data/", name_dataset, "/", sep = "")
path_images <- paste("images/", name_dataset, "/", sep = "")
mkdir(c(path_images))
path_images <- paste("images/", name_dataset, "/NSC_optimization/", sep = "")
mkdir(c(path_images))
path_results <- paste("results/", name_dataset, "/", sep = "")
mkdir(c(path_results))

## figures dimensions
figure_width <- 8
figure_height <- 8




#=============================================================================#




## ||||||||||||||||||||
# Load in Data ----
## ||||||||||||||||||||

load(paste(path_data, "mesh.RData", sep = ""))
load(paste(path_data, "analyzed_data.RData", sep = ""))
load(paste(path_results, "fPCA.RData", sep = ""))




#=============================================================================#



### Redefine RMSE

RMSE <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}



#=============================================================================#




## ||||||||||||||||||||
# Decide Optimal Number of NSCs ----
## ||||||||||||||||||||

### Get residuals and residuals_norm
norm <- RMSE(counts, 0)
residuals <- fPCA_residuals(counts, scores, loadings_locs)
residuals_norm <- list()
residuals_norm[["N"]] <- 1
residuals_norm <- c(residuals_norm, residuals/norm)


### Plot
custom_theme <- theme(
  plot.title = element_text(size = 20, face = "bold"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14)
)

plot <- plot.nComp_selection(residuals_norm, 9, 6, theme_custom = custom_theme)
grid.draw(plot)


### Save
filepath <- paste(path_images, "Optimal_nComp_Selection.pdf", sep = "")
ggsave(filename = filepath, plot = plot, 
       width = figure_width*3, height = figure_height*2, 
       dpi = "print", units = "cm")





## Decide Optimal Number of NSCs ----
## ||||||||||||||||||||

### Save final nComp_opt
nComp_opt <- 4  
filepath <- file.path(path_results, "nComp_opt.RData")
save(nComp_opt, file=filepath)

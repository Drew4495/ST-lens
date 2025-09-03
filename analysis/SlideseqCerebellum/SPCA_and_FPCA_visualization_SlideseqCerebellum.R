# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Plotting Slideseq Cerebellum dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
PLOT <- TRUE
fdaPDE <- TRUE

wd <- getwd()
setwd(wd)
source("src/utils/cat.R")




#=============================================================================#




## ||||||||||||||||||||
# Load in Libraries and Functions ----
## ||||||||||||||||||||


## Libraries ----
## ||||||||||||||||||||

cat.section_title("Libraries")

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
library(EnhancedVolcano)

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
library(data.table)

#Stats
library(car)



## Functions ----
## ||||||||||||||||||||

cat.section_title("Functions")

source("src/plots.R")
source("src/meshing.R")
source("src/fdaPDE.R")
source("src/clustering.R")

source("src/src_alignment.R")
source("src/src_gene_stats_tests.R")

source("src/utils/geometry.R")





#=============================================================================#




## ||||||||||||||||||||
# Global variables ----
## ||||||||||||||||||||
cat.section_title("Global variables")

# Directories
name.dataset <- "SlideseqCerebellum"
directory.initial_data <- glue("data/{name.dataset}/")
directory.results <- glue("results/{name.dataset}/")
directory.images <- glue("images/{name.dataset}/test_images/")
directory.images.manuscript <- glue("images/{name.dataset}/manuscript/")


# Code flow control
RUN <- list()
RUN[["Truth-SPCA-fPCA"]] <- TRUE
RUN[["High-Res fPCA"]] <- TRUE
RUN[["ARI/CHAOS Scoring"]] <- TRUE #This takes a long time. Should only do one and save results
RUN[["ARI/CHAOS Vs. NSCs"]] <- TRUE
RUN[["Clustering Vs. NSCs"]] <- TRUE
RUN[["Spatial Comps"]] <- TRUE
RUN[["Gene by Cluster"]] <- TRUE

# Run Manuscript Figures?
MAN_RUN <- TRUE

# Set color map
qual_cmap <- "Accent"
div_cmap <- coolwarm(200)
accent_colors <- brewer.pal(8, qual_cmap) #8 colors in Accent
accent_order <- c(7,5,1,4,3,6,8,2)
qual_custom_colors_gt <- accent_colors[accent_order]
qual_custom_colors_gt <- c(qual_custom_colors_gt, "#1b9e77")
accent_order <- c(5,3,9,6,4,7,8,1,2)
qual_custom_colors_gt <- qual_custom_colors_gt[accent_order]
qual_custom_colors_gt[9] <- "#83AEE3"




#=============================================================================#




## ||||||||||||||||||||
# Load Data ----
## ||||||||||||||||||||
# Original Data
load(glue("{directory.initial_data}data.RData"))

#SPCA results
load(glue("{directory.results}SPCA_results_{name.dataset}.RData"))

#fPCA results and data
load(glue("{directory.initial_data}mesh.RData"))
load(glue("{directory.initial_data}analyzed_data.RData"))
load(glue("{directory.results}fPCA.RData"))
load(glue("{directory.results}nComp_opt.RData"))
load(glue("{directory.results}fPCA_clustering_results_GCV_NoMean_6.RData"))
#load(glue("{directory.results}fPCA_clustering_results_GCV_NoMean_5.RData"))


#Defining PSEUDO true labels for fPCA with named atomic vector
pseudotruth_labels_fPCA <- cluster_labels_list_GCV_NoMean[["knearest_166"]]




#=============================================================================#




## ||||||||||||||||||||
# Alignment of Clustering Labels ----
## ||||||||||||||||||||

# fPCA LR alignment: already aligned

#SPCA alignment of 9 clusters based on Pseudotruth
names(SPCA_9clusterlabels) <- rownames(locations_final_SPCA)
SPCA_9clusters_aligned_list <- get_aligned_labels(clustering_labels = SPCA_9clusterlabels[names(pseudotruth_labels_fPCA)], 
                                                  ground_truth = pseudotruth_labels_fPCA, return_mapping=TRUE)
SPCA_gt_mapping <- SPCA_9clusters_aligned_list[["optimal_permutation"]]
aligned_SPCA_labels_9clusters <- align_with_mapping(SPCA_9clusterlabels, SPCA_gt_mapping)




#=============================================================================#




## ||||||||||||||||||||
# Scoring of Clustering Labels ----
## ||||||||||||||||||||

#fPCA - LR labels
fPCA_CHAOS <- fx_CHAOS(as.character(pseudotruth_labels_fPCA), locations)




#=============================================================================#




## ||||||||||||||||||||
# PLOT: Truth-SPCA-fPCA ----
## ||||||||||||||||||||

cat.subsection_title("Truth-SPCA-fPCA")

tic()
if(RUN[["Truth-SPCA-fPCA"]]){
  
  #Set some plotting variables
  scale_factor <- 5.5
  width = 0.86*scale_factor
  height = 0.86*scale_factor
  pt_size = 0.1
  tit_size = 6*scale_factor
  ax_size = 4*scale_factor
  
  #SPCA
  num_labels_SPCA <- length(unique(aligned_SPCA_labels_9clusters))
  qual_custom_colors_SPCA <- qual_custom_colors_gt[1:num_labels_SPCA]
  plot_SPCA_man <- plot.points(locations=locations_final_SPCA, data=aligned_SPCA_labels_9clusters, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_SPCA)
  plot_SPCA_man = plot_SPCA_man +
    ggtitle("SpatialPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(SPCA_CHAOS_9,3)}"))
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_SPCA.pdf"), 
         plot=plot_SPCA_man, width=width, height=height, units="in")
  print(plot_SPCA_man) 
  
  #fPCA - 9 clusters
  num_labels_FPCA <- length(unique(pseudotruth_labels_fPCA))
  qual_custom_colors_FPCA <- qual_custom_colors_gt[1:num_labels_FPCA]
  plot_FPCA_man <- plot.points(locations=locations, data=pseudotruth_labels_fPCA, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_FPCA)
  plot_FPCA_man = plot_FPCA_man +
    ggtitle("fPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(fPCA_CHAOS,3)}"))
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA.pdf"), 
         plot=plot_FPCA_man, width=width, height=height, units="in")
  print(plot_FPCA_man) 
  
  
}





#=============================================================================#




## ||||||||||||||||||||
# High-Res fPCA ----
## ||||||||||||||||||||

if(RUN[["High-Res fPCA"]]){
  
  #Set HR labels
  cluster_labels_HR.fPCA <- cluster_labels_HR_list_GCV_NoMean[["knearest_166"]]
  
  ## clustering grid
  grid_step <- 25
  grid <- square_grid(SpatialPoints(locations)@bbox, grid_step, seed_point = seed_point)
  grid <- grid[!is.na(over(SpatialPoints(grid), lattice$domain)),]
  
  #Set some plotting variables
  scale_factor <- 2.5
  width = 2.4*scale_factor
  height = 2.4*scale_factor
  pt_size = 0.1
  tit_size = 16.74*scale_factor
  ax_size = 11.16*scale_factor
  
  #fPCA 9 clusters
  #qual_custom_colors_FPCA <- qual_custom_colors_gt[c(3,2,5,4,1,9,6,8,7)]
  num_labels_FPCA <- length(unique(pseudotruth_labels_fPCA))
  qual_custom_colors_FPCA <- qual_custom_colors_gt[1:num_labels_FPCA]
  plot_FPCA.HR_man <- plot.points(locations=grid, data=cluster_labels_HR.fPCA, size=pt_size, discrete=TRUE, 
                                  custom_colors = qual_custom_colors_FPCA)
  plot_FPCA.HR_man = plot_FPCA.HR_man +
    ggtitle("fPCA - High Res") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab("")
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA_HR.pdf"), 
         plot=plot_FPCA.HR_man, width=width, height=height, units="in")
  print(plot_FPCA.HR_man) 
  
}

#=============================================================================#




## ||||||||||||||||||||
# PLOT: Spatial Comps ----
## ||||||||||||||||||||

cat.subsection_title("Spatial Comps")

#Plot comps for manuscript
if (MAN_RUN){
  #Select comps
  num_clusters <- 9
  selected_comps <- seq(num_clusters)
  ncol = num_clusters
  FPCA_comp_colors <- coolwarm(200)
  SPCA_comp_colors <- rev(FPCA_comp_colors)
  scale_factor <- 1.2
  pt_size = 0.1
  height = 15/4 * scale_factor
  width = height*ncol
  tit_size = 16.74*scale_factor
  ax_size = 11.16*scale_factor
  
  #Create formatting
  manuscript_Scomp_plot_settings <- function(){
    manuscript_plot_settings <- theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(hjust=0.5, size=30, family="sans", face='bold'),
      )
  }
  
  #FPCA HR selected comps without scaling
  plot_Scomps <- plot.Scomps(locations=grid, loadings=loadings_HR, 
                             size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix="", is_equal_scale=FALSE, ncol=ncol, 
                             colormap_continuous = FPCA_comp_colors)
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA_SpatiaComps_Selected_Unscaled.pdf"),
         plot=plot_Scomps, width=width, height=height, units='in')
  grid.arrange(plot_Scomps)
  
  #FPCA HR selected comps with scaling
  plot_Scomps_scaled <- plot.Scomps(locations=grid, loadings=loadings_HR, 
                                    size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                    colormap_continuous = FPCA_comp_colors)
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA_SpatiaCompsSelected_Scaled.pdf"),
         plot=plot_Scomps_scaled, width=width, height=height, units='in')
  grid.arrange(plot_Scomps_scaled)
  
  #format SPCA loadings to match locations
  SPCA_loadings <- SPCA_obj@SpatialPCs
  
  #SPCA selected comps without scaling
  plot_Scomps_SPCA <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_loadings)[,1:9],
                                  size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                  title_prefix="", is_equal_scale=FALSE, ncol=ncol,
                                  colormap_continuous = SPCA_comp_colors)
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_SPCA_SpatiaCompsSelected_Unscaled.pdf"),
         plot=plot_Scomps_SPCA, width=width, height=height, units='in')
  grid.arrange(plot_Scomps_SPCA)
  
  #SPCA selected comps with scaling
  plot_Scomps_SPCA_scaled <- plot.Scomps2(locations=locations_final_SPCA, loadings=t(SPCA_loadings)[ ,1:9],
                                          size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                          title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                          colormap_continuous = SPCA_comp_colors)
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_SPCA_SpatiaCompsSelected_Scaled.pdf"),
         plot=plot_Scomps_SPCA_scaled, width=width, height=height, units='in')
  grid.arrange(plot_Scomps_SPCA_scaled)
  
}






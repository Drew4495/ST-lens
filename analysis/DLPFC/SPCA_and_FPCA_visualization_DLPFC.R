# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Plotting DLPFC dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
PLOT <- TRUE
fdaPDE <- TRUE

setwd(getwd())
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
directory.initial_data <- "data/DLPFC/"
directory.results <- "results/DLPFC/"
directory.images <- "images/DLPFC/test_images/"
directory.images.manuscript <- "images/DLPFC/manuscript/"
name.dataset <- "DLPFC_sample9"

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
accent_order <- c(7,5,1,3,4,6,8,2)
qual_custom_colors_gt <- accent_colors[accent_order]




#=============================================================================#




## ||||||||||||||||||||
# Load Data ----
## ||||||||||||||||||||
# Original Data
load(glue("{directory.initial_data}data.RData"))

#SPCA results
#load(glue("{directory.results}SPCA_results_{name.dataset}.RData"))
load(glue("{directory.results}SPCA_allresults_NEW_{name.dataset}.RData"))

#fPCA results and data
load(glue("{directory.initial_data}analyzed_data.RData"))
load(glue("{directory.results}fPCA.RData"))
load(glue("{directory.results}nComp_opt.RData"))
load(glue("{directory.results}fPCA_clustering_results_GCV_NoMean_3.RData"))
load(glue("{directory.results}fPCA_clustering_results_ALIGNED_GCV_NoMean_3.RData"))

#Defining true labels for fPCA with named atomic vector
true_labels_named_vector <- true_labels$true_label
names(true_labels_named_vector) <- rownames(true_labels)
true_labels_fPCA <- true_labels_named_vector[rownames(locations)]

#Define labels of choice
fPCA_labels_aligned <- cluster_labels_list_aligned_GCV_NoMean[["knearest_240"]]




#=============================================================================#




## ||||||||||||||||||||
# Alignment of Clustering Labels ----
## ||||||||||||||||||||

# fPCA LR alignment: already aligned

# SPCA alignment: already aligned




#=============================================================================#




## ||||||||||||||||||||
# Scoring of Clustering Labels ----
## ||||||||||||||||||||
#fPCA - LR labels
fPCA_ARI <- adjustedRandIndex(fPCA_labels_aligned, true_labels_fPCA)
fPCA_CHAOS <- fx_CHAOS(as.character(true_labels_fPCA), locations)




#=============================================================================#




## ||||||||||||||||||||
# PLOT: Truth-SPCA-fPCA ----
## ||||||||||||||||||||

cat.subsection_title("Truth-SPCA-fPCA")

tic()
if(RUN[["Truth-SPCA-fPCA"]]){
  
  #Set some plotting variables
  width = 0.86*3
  height = 0.86*3
  pt_size = 0.4
  tit_size = 6*3
  ax_size = 4*3
  
  #Gound Truth
  plot_gt_man <- plot.points(locations=locations, data=true_labels, size=pt_size, discrete=TRUE, 
                             custom_colors=qual_custom_colors_gt)
  plot_gt_man = plot_gt_man +
    ggtitle("Ground Truth") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5))
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_GroundTruth.pdf"), 
         plot=plot_gt_man, width=width, height=height, units="in")
  print(plot_gt_man) 
  
  #SPCA
  num_labels_SPCA <- length(unique(aligned_SPCA_labels))
  qual_custom_colors_SPCA <- qual_custom_colors_gt[1:num_labels_SPCA]
  plot_SPCA_man <- plot.points(locations=locations_final_SPCA, data=aligned_SPCA_labels, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_SPCA)
  plot_SPCA_man = plot_SPCA_man +
    ggtitle("SpatialPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("ARI: {round(SPCA_ARI,2)}  |  CHAOS: {round(SPCA_CHAOS,2)}"))
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_SPCA.pdf"), 
         plot=plot_SPCA_man, width=width, height=height, units="in")
  print(plot_SPCA_man) 
  
  #fPCA
  num_labels_FPCA <- length(unique(fPCA_labels_aligned))
  qual_custom_colors_FPCA <- qual_custom_colors_gt[1:num_labels_FPCA]
  plot_FPCA_man <- plot.points(locations=locations, data=fPCA_labels_aligned, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_FPCA)
  plot_FPCA_man = plot_FPCA_man +
    ggtitle("fPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("ARI: {round(fPCA_ARI,2)}  |  CHAOS: {round(fPCA_CHAOS,2)}"))
  ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA.pdf"), 
         plot=plot_FPCA_man, width=width, height=height, units="in")
  print(plot_FPCA_man) 
  
  
}




#=============================================================================#



## ||||||||||||||||||||
# PLOT: Spatial Comps ----
## ||||||||||||||||||||

cat.subsection_title("Spatial Comps")

tic()
if(RUN[["Spatial Comps"]]){
  
  #Create formatting
  manuscript_Scomp_plot_settings <- function(){
    manuscript_plot_settings <- theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(hjust=0.5, size=16, family="sans", face='bold'),
      )
  }
  
  #Plot comps without scaling
  plot_Scomps <- plot.Scomps(locations=locations, loadings=loadings_locs, 
                             size=0.25, extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix="", is_equal_scale=FALSE)
  grid.arrange(plot_Scomps)
  ggsave(filename=glue("{directory.images.manuscript}HER2_FPCA_SpatialComps.pdf"), 
         plot=plot_Scomps, width=4, height=4, units="in")
  
  #Plot comps with scaling
  plot_Scomps_scaled <- plot.Scomps(locations=locations, loadings=loadings_locs, 
                                    size=0.25, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=TRUE)
  grid.arrange(plot_Scomps_scaled)
  ggsave(filename=glue("{directory.images.manuscript}HER2_FPCA_SpatialComps_Scaled.pdf"), 
         plot=plot_Scomps_scaled, width=4, height=4, units="in")
  
  #Plot comps for manuscript
  if (MAN_RUN){
    #Select comps
    selected_comps <- seq(7)
    ncol = 7
    FPCA_comp_colors <- coolwarm(200)
    SPCA_comp_colors <- rev(FPCA_comp_colors)
    pt_size = 0.2
    height = 2
    width = 2*ncol
    
    
    #FPCA selected comps without scaling
    plot_Scomps <- plot.Scomps(locations=locations, loadings=loadings_locs[ ,selected_comps], 
                               size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                               title_prefix="", is_equal_scale=FALSE, ncol=ncol, 
                               colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA_SpatiaComps_Selected_Unscaled.pdf"),
           plot=plot_Scomps, width=width, height=height, units='in')
    grid.arrange(plot_Scomps)
    
    #FPCA selected comps with scaling
    plot_Scomps_scaled <- plot.Scomps(locations=locations, loadings=loadings_locs[ ,selected_comps], 
                                      size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                      title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                      colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_FPCA_SpatiaCompsSelected_Scaled.pdf"),
           plot=plot_Scomps_scaled, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_scaled)
    
    #format SPCA loadings to match locations
    SPCA_loadings <- SPCA_obj@SpatialPCs
    SPCA_loadings <- SPCA_loadings[ ,rownames(locations_final_SPCA)]
    
    #SPCA selected comps without scaling
    plot_Scomps_SPCA <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_loadings)[ ,selected_comps],
                                    size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=FALSE, ncol=ncol,
                                    colormap_continuous = SPCA_comp_colors)
    ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_SPCA_SpatiaCompsSelected_Unscaled.pdf"),
           plot=plot_Scomps_SPCA, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_SPCA)
    
    #SPCA selected comps with scaling
    plot_Scomps_SPCA_scaled <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_loadings)[ ,selected_comps],
                                           size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                           title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                           colormap_continuous = SPCA_comp_colors)
    ggsave(filename=glue("{directory.images.manuscript}{name.dataset}_SPCA_SpatiaCompsSelected_Scaled.pdf"),
           plot=plot_Scomps_SPCA_scaled, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_SPCA_scaled)
    
  }
  
}






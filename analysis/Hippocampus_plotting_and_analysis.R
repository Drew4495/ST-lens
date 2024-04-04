# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Plotting HER2 dataset %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

VERBOSE <- TRUE
PLOT <- TRUE
fdaPDE <- FALSE

wd <- "/Users/aburns/Codebook/Projects/spatial_transcriptomics_GRSVD/Paper/ST-lens"
setwd(wd)
source("src/cat_utilities.R")

cat.script_title("Plotting Hippocampus dataset")



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
source("src/fdaPDE.R")
source("src/clustering.R")

source("src/src_alignment.R")



# |||||||||||||||||||||
# Global variables ----
# |||||||||||||||||||||

cat.section_title("Global variables")

# Directories
directory.initial_data <- "data/SlideseqV2Hippocampus/"
directory.results <- "results/SlideseqV2Hippocampus/"
directory.images <- "images/SlideseqV2Hippocampus/"
directory.images.manuscript <- "images/SlideseqV2Hippocampus/manuscript/"
name.dataset <- "SlideseqV2Hippocampus"

# Code flow control
RUN <- list()
RUN[["Truth-SPCA-fPCA"]] <- TRUE
RUN[["High-Res fPCA"]] <- TRUE
RUN[["ARI/CHAOS Scoring"]] <- FALSE #This takes a long time. Should only do one and save results
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
qual_custom_colors_gt <- c(qual_custom_colors_gt, "#1b9e77")




# |||||||||||||||||||||
# Data ----
# |||||||||||||||||||||

# Original Data
load(glue("{directory.initial_data}{name.dataset}.RData"))

#SPCA results
#load(glue("{directory.results}SPCA_results_{name.dataset}.RData"))

#fPCA results
load(glue("{directory.results}{name.dataset}_processed.RData"))
load(glue("{directory.results}{name.dataset}_mesh.RData"))
load(glue("{directory.results}{name.dataset}_mean.RData"))
load(glue("{directory.results}{name.dataset}_nComp_selection.RData"))
#load(glue("{directory.results}{name.dataset}_components.RData")) #No such file or directory
load(glue("{directory.results}{name.dataset}_clusters_on_HR_grid.RData"))
load(glue("{directory.results}grid.RData")) #locations for HR dataset

#Defining PSEUDO true labels for fPCA with named atomic vector
true_labels_fPCA <- cluster_labels_downsamples #Pseudotruth is just LR fPCA labels



# |||||||||||||||||||||
# Alignment of Clustering Labels ----
# |||||||||||||||||||||

# fPCA LR alignment
fPCA_labels_aligned <- cluster_labels_downsamples #already aligned since its pseudotruth

#fPCA HR alignment
fPCA.HR_labels_aligned <- cluster_labels_HR



# |||||||||||||||||||||
# Scoring of Clustering Labels ----
# |||||||||||||||||||||
#fPCA - LR labels
fPCA_CHAOS <- fx_CHAOS(as.character(true_labels_fPCA), locations.final)
fPCA_HR_CHAOS <- fx_CHAOS(as.character(fPCA.HR_labels_aligned), grid)



# |||||||||||||||||||||
# Truth-SPCA-fPCA
# |||||||||||||||||||||

cat.subsection_title("Truth-SPCA-fPCA")

tic()
if(RUN[["Truth-SPCA-fPCA"]]){
  
  #Set some plotting variables
  scale_factor <- 9
  width = 0.86*scale_factor
  height = 0.86*scale_factor
  pt_size = 0.1
  tit_size = 6*scale_factor
  ax_size = 4*scale_factor
  
  #fPCA
  #qual_custom_colors_FPCA <- qual_custom_colors_gt[c(3,2,5,4,1,9,6,8,7)]
  plot_FPCA_man <- plot.points(locations=locations.final, data=fPCA_labels_aligned, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_gt)
  plot_FPCA_man = plot_FPCA_man +
    ggtitle("fPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(fPCA_CHAOS,4)}"))
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA.pdf"), 
         plot=plot_FPCA_man, width=width, height=height, units="in")
  print(plot_FPCA_man) 
  
  
}



# |||||||||||||||||||||
# High-Res fPCA ----
# |||||||||||||||||||||

if(RUN[["High-Res fPCA"]]){
  
  #Set some plotting variables
  scale_factor <- 1.5
  width = 2.4*scale_factor
  height = 2.4*scale_factor
  pt_size = 0.3
  tit_size = 16.74*scale_factor
  ax_size = 11.16*scale_factor
  
  #fPCA - HR
  #qual_custom_colors_FPCA <- qual_custom_colors_gt[c(3,2,5,4,1,9,6,8,7)]
  plot_FPCA.HR_man <- plot.points(locations=grid, data=fPCA.HR_labels_aligned, size=pt_size, discrete=TRUE, 
                                  custom_colors = qual_custom_colors_gt)
  plot_FPCA.HR_man = plot_FPCA.HR_man +
    ggtitle("fPCA - High Res") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(fPCA_HR_CHAOS,4)}"))
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_HR.pdf"), 
         plot=plot_FPCA.HR_man, width=width, height=height, units="in")
  print(plot_FPCA.HR_man) 
  
  
  
  
}



# |||||||||||||||||||||
# Spatial Comps ----
# |||||||||||||||||||||

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
  
  #Create HR loadings. Takes up less memory for plotting.
  loadings_HR <- NULL
  for(h in 1:ncol(loadings)){
    loadings_HR <- cbind(loadings_HR, field.eval(grid, loadings_nodes[,h], mesh))
  }
  
  #plot settings - LR
  pt_size = 0.1
  height = 25
  width = 25
  
  #Plot ALL comps without scaling
  plot_Scomps <- plot.Scomps(locations=locations.final, loadings=loadings, 
                             size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix="", is_equal_scale=FALSE)
  grid.arrange(plot_Scomps)
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatialComps.pdf"), 
         plot=plot_Scomps, width=width, height=height, units="in")
  
  #Plot ALL comps with scaling
  plot_Scomps_scaled <- plot.Scomps(locations=locations.final, loadings=loadings, 
                                    size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=TRUE)
  grid.arrange(plot_Scomps_scaled)
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatialComps_Scaled.pdf"), 
         plot=plot_Scomps_scaled, width=width, height=height, units="in")
  
  #plot settings - HR
  pt_size = 0.3
  height = 15
  width = 15
  
  #Plot ALL comps without scaling - HR
  plot_Scomps <- plot.Scomps(locations=grid, loadings=loadings_HR, 
                             size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix="", is_equal_scale=FALSE)
  grid.arrange(plot_Scomps)
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatialComps_HR.pdf"), 
         plot=plot_Scomps, width=width, height=height, units="in")
  
  #Plot ALL comps with scaling - HR
  plot_Scomps_scaled <- plot.Scomps(locations=grid, loadings=loadings_HR, 
                                    size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=TRUE)
  grid.arrange(plot_Scomps_scaled)
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatialComps_Scaled_HR.pdf"), 
         plot=plot_Scomps_scaled, width=width, height=height, units="in")
  
  #Plot comps for manuscript
  if (MAN_RUN){
    #Select comps
    num_clusters <- 8
    selected_comps <- seq(num_clusters)
    ncol = num_clusters
    FPCA_comp_colors <- coolwarm(200)
    SPCA_comp_colors <- rev(FPCA_comp_colors)
    pt_size = 0.15
    height = 7
    width = height*ncol
    
    #Create formatting
    manuscript_Scomp_plot_settings <- function(){
      manuscript_plot_settings <- theme_void() +
        theme(legend.position="none", 
              plot.title=element_text(hjust=0.5, size=60, family="sans", face='bold'),
        )
    }
    
    #FPCA selected comps without scaling
    plot_Scomps <- plot.Scomps(locations=locations.final, loadings=loadings[ ,selected_comps], 
                               size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                               title_prefix="", is_equal_scale=FALSE, ncol=ncol, 
                               colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatiaComps_Selected_Unscaled.pdf"),
           plot=plot_Scomps, width=width, height=height, units='in', limitsize=FALSE)
    grid.arrange(plot_Scomps)
    
    #FPCA selected comps with scaling
    plot_Scomps_scaled <- plot.Scomps(locations=locations.final, loadings=loadings[ ,selected_comps], 
                                      size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                      title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                      colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatiaCompsSelected_Scaled.pdf"),
           plot=plot_Scomps_scaled, width=width, height=height, units='in', limitsize=FALSE)
    grid.arrange(plot_Scomps_scaled)
    
    #Select comps
    num_clusters <- 8
    selected_comps <- seq(num_clusters)
    ncol = num_clusters
    FPCA_comp_colors <- coolwarm(200)
    SPCA_comp_colors <- rev(FPCA_comp_colors)
    pt_size = 0.3
    height = 15/4
    width = height*ncol
    
    #Create formatting
    manuscript_Scomp_plot_settings <- function(){
      manuscript_plot_settings <- theme_void() +
        theme(legend.position="none", 
              plot.title=element_text(hjust=0.5, size=30, family="sans", face='bold'),
        )
    }
    
    #FPCA selected comps without scaling - HR
    plot_Scomps <- plot.Scomps(locations=grid, loadings=loadings_HR[ ,selected_comps], 
                               size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                               title_prefix="", is_equal_scale=FALSE, ncol=ncol, 
                               colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatiaComps_Selected_Unscaled_HR.pdf"),
           plot=plot_Scomps, width=width, height=height, units='in')
    grid.arrange(plot_Scomps)
    
    #FPCA selected comps with scaling - HR
    plot_Scomps_scaled <- plot.Scomps(locations=grid, loadings=loadings_HR[ ,selected_comps], 
                                      size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                      title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                      colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatiaCompsSelected_Scaled_HR.pdf"),
           plot=plot_Scomps_scaled, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_scaled)
    
  }
  
}


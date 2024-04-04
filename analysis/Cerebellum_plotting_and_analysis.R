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

cat.script_title("Plotting Cerebellum dataset")



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
directory.initial_data <- "data/SlideseqCerebellum/"
directory.results <- "results/slideseqCerebellum/"
directory.images <- "images/SlideseqCerebellum/"
directory.images.manuscript <- "images/SlideseqCerebellum/manuscript/"
name.dataset <- "SlideseqCerebellum"

# Code flow control
RUN <- list()
RUN[["fPCA: 8 clusters"]] <- FALSE
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
load(glue("{directory.results}SPCA_results_{name.dataset}.RData"))

#fPCA results
load(glue("{directory.results}{name.dataset}_processed.RData"))
load(glue("{directory.results}{name.dataset}_mesh.RData"))
load(glue("{directory.results}{name.dataset}_mean.RData"))
load(glue("{directory.results}{name.dataset}_nComp_selection.RData"))
#load(glue("{directory.results}{name.dataset}_components.RData")) #No such file or directory
load(glue("{directory.results}{name.dataset}_clusters_on_HR_grid.RData"))

#Defining PSEUDO true labels for fPCA with named atomic vector
true_labels_fPCA <- cluster_labels_downsamples #Pseudotruth is just LR fPCA labels



# |||||||||||||||||||||
# Alignment of Clustering Labels ----
# |||||||||||||||||||||

# fPCA LR alignment
fPCA_labels_aligned <- cluster_labels_downsamples #already aligned since its pseudotruth

#fPCA HR alignment
fPCA.HR_labels_aligned <- cluster_labels_refined_HR

#SPCA alignment of 8 clusters based on Pseudotruth
#names(aligned_SPCA_labels) <- rownames(locations_final_SPCA)
#SPCA_aligned_list <- get_aligned_labels(clustering_labels = aligned_SPCA_labels[names(fPCA_labels_aligned)], 
#                                         ground_truth = fPCA_labels_aligned, return_mapping=TRUE)
#SPCA_gt_mapping <- SPCA_aligned_list[["optimal_permutation"]]
#aligned_SPCA_labels <- align_with_mapping(aligned_SPCA_labels, SPCA_gt_mapping)


#SPCA alignment of 9 clusters based on Pseudotruth
names(SPCA_9clusterlabels) <- rownames(locations_final_SPCA)
SPCA_9clusters_aligned_list <- get_aligned_labels(clustering_labels = SPCA_9clusterlabels[names(fPCA_labels_aligned)], 
                                          ground_truth = fPCA_labels_aligned, return_mapping=TRUE)
SPCA_gt_mapping <- SPCA_9clusters_aligned_list[["optimal_permutation"]]
aligned_SPCA_labels_9clusters <- align_with_mapping(SPCA_9clusterlabels, SPCA_gt_mapping)


# |||||||||||||||||||||
# fPCA: 8 clusters ----
# |||||||||||||||||||||

if(RUN[["fPCA: 8 clusters"]]){
  #fPCA - LR
  clusternum = 8
  knearest <- round(sqrt(dim(SPCA_obj@SpatialPCs)[2]))
  fPCA_labels_8clusters <- louvain_clustering(clusternum, 
                                              latent_dat=t(loadings[ ,1:nComp_opt]),
                                              knearest = knearest)
  filepath <- glue("{directory.results}{name.dataset}_fPCA_labels_8clusters_LR.RData")
  save(fPCA_labels_8clusters, file=filepath)
  
  #fPCA - HR
  clusternum <- 8
  knearest <- round(sqrt(nrow(grid)))
  
  data.clustering <- NULL
  for(i in 1:nComp_opt){
    loading_HR <- field.eval(grid, loadings_nodes[,i], mesh)
    data.clustering <- cbind(data.clustering, loading_HR)
  }
  
  fPCA.HR_labels_8clusters <- louvain_clustering(clusternum = clusternum,
                                                 latent_dat = t(data.clustering[,1:nComp_opt]),
                                                 knearest = knearest)
  filepath <- glue("{directory.results}{name.dataset}_fPCA_labels_8clusters_HR.RData")
  save(fPCA_labels_8clusters, file=filepath)
  
}

#Load in above clustering.
filepath <- glue("{directory.results}{name.dataset}_fPCA_labels_8clusters_LR.RData")
load(filepath)
filepath <- glue("{directory.results}{name.dataset}_fPCA_labels_8clusters_HR.RData")
load(filepath)




# |||||||||||||||||||||
# Scoring of Clustering Labels ----
# |||||||||||||||||||||
### 9 clusters
#CHAOS
fPCA_CHAOS <- fx_CHAOS(as.character(true_labels_fPCA), locations.final)
fPCA_HR_CHAOS <- fx_CHAOS(as.character(fPCA.HR_labels_aligned), grid)
SPCA_CHAOS_9clusters <- fx_CHAOS(as.character(aligned_SPCA_labels_9clusters), locations_final_SPCA)

#LISI
fPCA_lisi <- compute_lisi(locations.final, data.frame("fPCA"=true_labels_fPCA), c('fPCA'))
fPCA_HR_lisi <- compute_lisi(grid, data.frame("fPCA"=fPCA.HR_labels_aligned), c('fPCA'))
SPCA_lisi <- compute_lisi(locations_final_SPCA, data.frame("SPCA"=aligned_SPCA_labels_9clusters), c('SPCA'))

#PAS
fPCA_PAS <- fx_PAS(true_labels_fPCA, locations.final)
fPCA_PAS_HR <- fx_PAS(fPCA.HR_labels_aligned, grid)
SPCA_PAS <- fx_PAS(aligned_SPCA_labels_9clusters, locations_final_SPCA)


### 8 clusters
#CHAOS
fPCA_CHAOS_8clusters <- fx_CHAOS(as.character(fPCA_labels_8clusters), locations.final)
#fPCA_HR_CHAOS_8clusters <- fx_CHAOS(as.character()
SPCA_CHAOS_8clusters <- fx_CHAOS(as.character(aligned_SPCA_labels), locations_final_SPCA)

#LISI
fPCA_lisi_8clusters <- compute_lisi(locations.final, data.frame("fPCA"=fPCA_labels_8clusters), c('fPCA'))
#fPCA_HR_lisi <- compute_lisi()
SPCA_lisi_8clusters <- compute_lisi(locations_final_SPCA, data.frame("SPCA"=aligned_SPCA_labels), c("SPCA"))

#PAS
fPCA_PAS_8clusters <- fx_PAS(fPCA_labels_8clusters, locations.final)
#fPCA_HR_PAS_8clusters <- fx_PAS()
SPCA_PAS_8clusters <-fx_PAS(aligned_SPCA_labels, locations_final_SPCA)


###Saving all scores because some take a lot of memory and time
filepath <- glue("{directory.initial_data}{name.dataset}_AllScores_FPCAandSPCA")
save(fPCA_CHAOS, fPCA_HR_CHAOS, SPCA_CHAOS_9clusters,
     fPCA_lisi, fPCA_HR_lisi, SPCA_lisi,
     fPCA_PAS, fPCA_PAS_HR, SPCA_PAS,
     fPCA_CHAOS_8clusters, SPCA_CHAOS_8clusters,
     fPCA_lisi_8clusters, SPCA_lisi_8clusters,
     fPCA_PAS_8clusters, SPCA_PAS_8clusters, file=filepath)



# |||||||||||||||||||||
# Truth-SPCA-fPCA
# |||||||||||||||||||||

cat.subsection_title("Truth-SPCA-fPCA")

tic()
if(RUN[["Truth-SPCA-fPCA"]]){
  
  #Set some plotting variables
  width = 0.86*6
  height = 0.86*6
  pt_size = 0.7
  tit_size = 6*6
  ax_size = 4*6

  
  #SPCA - 8 clusters
  plot_SPCA_man <- plot.points(locations=locations_final_SPCA, data=aligned_SPCA_labels, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_gt)
  plot_SPCA_man = plot_SPCA_man +
    ggtitle("SpatialPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(SPCA_CHAOS,4)}"))
  ggsave(filename=glue("{directory.images}{name.dataset}_SPCA_8clusters.pdf"), 
         plot=plot_SPCA_man, width=width, height=height, units="in")
  print(plot_SPCA_man) 
  
  #SPCA - 9 clusters
  qual_custom_colors_FPCA <- qual_custom_colors_gt[c(3,2,5,4,1,9,6,8,7)]
  num_labels_SPCA <- length(unique(aligned_SPCA_labels))
  qual_custom_colors_SPCA <- qual_custom_colors_gt[1:num_labels_SPCA]
  plot_SPCA_man <- plot.points(locations=locations_final_SPCA, data=aligned_SPCA_labels_9clusters, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_FPCA)
  plot_SPCA_man = plot_SPCA_man +
    ggtitle("SpatialPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(SPCA_CHAOS_9clusters,4)}"))
  ggsave(filename=glue("{directory.images}{name.dataset}_SPCA_9clusters.pdf"), 
         plot=plot_SPCA_man, width=width, height=height, units="in")
  print(plot_SPCA_man) 
  
  #fPCA - 8 clusters
  qual_custom_colors_FPCA_8clusters <- qual_custom_colors_gt[c(9,1,7,2,3,4,5,6,8)]
  plot_FPCA_man <- plot.points(locations=locations.final, data=fPCA_labels_8clusters, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_FPCA_8clusters)
  plot_FPCA_man = plot_FPCA_man +
    ggtitle("fPCA") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(fPCA_CHAOS_8clusters,4)}"))
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_8clusters.pdf"), 
         plot=plot_FPCA_man, width=width, height=height, units="in")
  print(plot_FPCA_man) 
  
  #fPCA - 9 clusters
  qual_custom_colors_FPCA <- qual_custom_colors_gt[c(3,2,5,4,1,9,6,8,7)]
  plot_FPCA_man <- plot.points(locations=locations.final, data=fPCA_labels_aligned, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_FPCA)
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
  scale_factor <- 1.45
  width = 2.4*scale_factor
  height = 2.4*scale_factor
  pt_size = 0.2
  tit_size = 16.74*scale_factor
  ax_size = 11.16*scale_factor
  
  #fPCA 8 clusters
  qual_custom_colors_FPCA_8clusters <- qual_custom_colors_gt[c(3,2,5,4,9,6,1,7,8)]
  plot_FPCA.HR_man <- plot.points(locations=grid, data=fPCA.HR_labels_8clusters, size=pt_size, discrete=TRUE, 
                                  custom_colors = qual_custom_colors_FPCA_8clusters)
  plot_FPCA.HR_man = plot_FPCA.HR_man +
    ggtitle("fPCA - High Res") +
    theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
          axis.title.x=element_text(size=ax_size, family="sans")) +
    xlab(glue("CHAOS: {round(fPCA_HR_CHAOS_8clusters,4)}"))
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_HR_8clusters.pdf"), 
         plot=plot_FPCA.HR_man, width=width, height=height, units="in")
  print(plot_FPCA.HR_man) 
  
  #fPCA 9 clusters
  qual_custom_colors_FPCA <- qual_custom_colors_gt[c(3,2,5,4,1,9,6,8,7)]
  plot_FPCA.HR_man <- plot.points(locations=grid, data=fPCA.HR_labels_aligned, size=pt_size, discrete=TRUE, 
                                  custom_colors = qual_custom_colors_FPCA)
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
# ARI/CHAOS Scoring ----
# |||||||||||||||||||||

## ||||||||||||||||||||
## NOTE:  I may want to put this in another script ----
## |||||||||||||||||||||

cat.subsection_title("ARI/CHAOS Scoring")

tic()
if(RUN[["ARI/CHAOS Scoring"]]){
  
  #initialize variables
  comp_seq <- seq(dim(loadings)[2])
  score_by_NSCs_LR_df <- data.frame(matrix(NA, nrow=length(comp_seq), ncol=3))
  score_by_NSCs_HR_df <- data.frame(matrix(NA, nrow=length(comp_seq), ncol=3))
  colnames(score_by_NSCs_LR_df) <- c("ARI", "CHAOS", "PAS")
  rownames(score_by_NSCs_LR_df) <- comp_seq
  colnames(score_by_NSCs_HR_df) <- c("ARI", "CHAOS", "PAS")
  rownames(score_by_NSCs_HR_df) <- comp_seq
  clusternum <- 6
  
  data.clustering <- NULL
  for(h in 1:length(comp_seq)){ 
    loading_HR <- field.eval(grid, loadings_nodes[,h], mesh)
    data.clustering <- cbind(data.clustering, loading_HR)
  }
  
  cluster_labels_LR_all <- as.data.frame(matrix(NA, nrow=dim(loadings)[1], ncol=length(comp_seq)))
  cluster_labels_HR_all <- as.data.frame(matrix(NA, nrow=dim(data.clustering)[1], ncol=length(comp_seq)))
  #loop through cluster combinations
  for (i in comp_seq){ 
    #Perform clustering for LR
    knearest <- round(sqrt(nrow(locations.final)))
    cluster_labels <- walktrap_clustering(clusternum = clusternum,
                                          latent_dat = t(loadings[ ,1:i]),
                                          knearest = knearest)
    names(cluster_labels) <- rownames(locations.final)
    cluster_labels_LR_all[ ,i] <- cluster_labels
    
    #select labels that don't correspond to undetermined ground truth
    names.locations_unknown <- names(true_labels_fPCA[true_labels_fPCA == 7])
    cluster_labels_no7 <- cluster_labels[!names(cluster_labels) %in% names.locations_unknown]
    
    #scoring for LR
    ARI <- adjustedRandIndex(true_labels_no7, cluster_labels_no7)
    CHAOS <- fx_CHAOS(cluster_labels_no7, locations.final)
    #PAS <- fx_PAS(cluster_labels_no7, locations.final)
    vals <- c(ARI, CHAOS, NA)
    score_by_NSCs_LR_df[i, ] <- vals
    cat("\n","LR", i, "done")
    
    
    # Clustering for HR
    knearest <- 200
    cluster_labels_HR <- walktrap_clustering(clusternum = clusternum,
                                             latent_dat = t(data.clustering[ ,1:i]),
                                             knearest = knearest)
    
    cluster_labels_HR_all[ ,i] <- cluster_labels_HR
    cluster_labels <- c()
    names.locations <- rownames(locations.final)
    len.HR <- length(cluster_labels_HR)
    for(name.l in names.locations){
      
      distances <- dist_point_from_points(locations.final[name.l,], grid)
      index.closest_point <- which.min(distances)
      cluster_labels[name.l] <- cluster_labels_HR[index.closest_point]
      
    }
    
    #select labels that don't correspond to undetermined ground truth
    names.locations_unknown <- names(true_labels_fPCA[true_labels_fPCA == 7])
    cluster_labels_no7 <- cluster_labels[!names(cluster_labels) %in% names.locations_unknown]
    
    #scoring for HR
    ARI <- adjustedRandIndex(true_labels_no7, cluster_labels_no7)
    CHAOS <- fx_CHAOS(cluster_labels_no7, locations.final)
    #PAS <- fx_PAS(cluster_labels_no7, locations.final)
    vals <- c(ARI, CHAOS, NA)
    score_by_NSCs_HR_df[i, ] <- vals
    cat("\n","HR", i, "done")
    
  }
  
  
  #Save score dataframes
  #fPCA_HR_ncomp3_ARI <- score_by_NSCs_HR_df$ARI[3]
  #fPCA_HR_ncomp3_CHAOS <- score_by_NSCs_HR_df$CHAOS[3]
  filepath <- glue("{directory.results}HER2_HRandLR_clustering_scores.RData")
  save(score_by_NSCs_HR_df, score_by_NSCs_LR_df, cluster_labels_LR_all,
       cluster_labels_HR_all, file=filepath)
  
  
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
            plot.title=element_text(hjust=0.5, size=30, family="sans", face='bold'),
      )
  }
  
  #Create HR loadings. Takes up less memory for plotting.
  loadings_HR <- NULL
  for(h in 1:ncol(loadings)){
    loadings_HR <- cbind(loadings_HR, field.eval(grid, loadings_nodes[,h], mesh))
  }
  
  #plot settings - LR
  pt_size = 0.1
  height = 15
  width = 15
  
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
    num_clusters <- 9
    selected_comps <- seq(num_clusters)
    ncol = num_clusters
    FPCA_comp_colors <- coolwarm(200)
    SPCA_comp_colors <- rev(FPCA_comp_colors)
    pt_size = 0.1
    height = 15/4
    width = height*ncol
    
    #Create formatting
    manuscript_Scomp_plot_settings <- function(){
      manuscript_plot_settings <- theme_void() +
        theme(legend.position="none", 
              plot.title=element_text(hjust=0.5, size=30, family="sans", face='bold'),
        )
    }
    
    #FPCA selected comps without scaling
    plot_Scomps <- plot.Scomps(locations=locations.final, loadings=loadings[ ,selected_comps], 
                               size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                               title_prefix="", is_equal_scale=FALSE, ncol=ncol, 
                               colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatiaComps_Selected_Unscaled.pdf"),
           plot=plot_Scomps, width=width, height=height, units='in')
    grid.arrange(plot_Scomps)
    
    #FPCA selected comps with scaling
    plot_Scomps_scaled <- plot.Scomps(locations=locations.final, loadings=loadings[ ,selected_comps], 
                                      size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                      title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                      colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatiaCompsSelected_Scaled.pdf"),
           plot=plot_Scomps_scaled, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_scaled)
    
    #format SPCA loadings to match locations
    SPCA_loadings <- SPCA_obj@SpatialPCs
    
    #SPCA selected comps without scaling
    plot_Scomps_SPCA <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_loadings)[ ,selected_comps],
                                    size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=FALSE, ncol=ncol,
                                    colormap_continuous = SPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_SPCA_SpatiaCompsSelected_Unscaled.pdf"),
           plot=plot_Scomps_SPCA, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_SPCA)
    
    #SPCA selected comps with scaling
    plot_Scomps_SPCA_scaled <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_loadings)[ ,selected_comps],
                                           size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                           title_prefix="", is_equal_scale=TRUE, ncol=ncol,
                                           colormap_continuous = SPCA_comp_colors)
    ggsave(filename=glue("{directory.images}{name.dataset}_SPCA_SpatiaCompsSelected_Scaled.pdf"),
           plot=plot_Scomps_SPCA_scaled, width=width, height=height, units='in')
    grid.arrange(plot_Scomps_SPCA_scaled)
    
    #Select comps
    num_clusters <- 9
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


# |||||||||||||||||||||
# Violin Plot of scores
# |||||||||||||||||||||
fPCA_lisi <- as.data.table(fPCA_lisi, keep.rownames="RowLabel")
SPCA_lisi <- as.data.table(SPCA_lisi, keep.rownames="RowLabel")
lisi_df <- merge(fPCA_lisi, SPCA_lisi, by="RowLabel", all=TRUE)
lisi_df <- melt(lisi_df[ ,2:3], measure_vars=c('fPCA', 'SPCA'), variable.name='Method', value.name='Lisi')
  
ggplot(lisi_df, aes(x=Method, y=Lisi)) +
  geom_violin(trim=FALSE) + 
  labs(x="Method", y="Lisi scores")



# |||||||||||||||||||||
# De-noised Gene Expression
# |||||||||||||||||||||

###Denoised gene expression

#Identify gene
top_gene_idx <- which(scores[ ,1] == max(scores[ ,1]))
top_gene_name <- rownames(counts.significant)[top_gene_idx]

#Plotting params
width = 0.86*6
height = 0.86*6
pt_size = 0.3
tit_size = 6*6
ax_size = 4*6
Reds <- c("#f5f5dc", "darkred")

#Plot gene expression after SCTransform function
plot.unnoised <- plot.points(locations=locations.final, counts.centered[top_gene_idx, ], 
            size=pt_size, cmap_continuous=Reds)
plot.unnoised = plot.unnoised +
  ggtitle("Preprocessed") +
  theme_void() +
  theme(legend.position="right", legend.title = element_blank(), 
        plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
        axis.title.x=element_text(size=ax_size, family="sans")) +
  xlab(glue(top_gene_name))

ggsave(filename=glue("{directory.images}{name.dataset}_{top_gene_name}_unnoised.pdf"), 
       plot=plot.unnoised, width=width, height=height, units="in")
print(plot.unnoised) 


#Plot denoised expression
denoised_counts.significant <- scores[ ,1:9] %*% t(loadings[ ,1:9])
plot.denoised_pcp4 <- plot.points(locations=locations.final, denoised_counts.significant[top_gene_idx, ], 
            size=pt_size, cmap_continuous=Reds)
plot.denoised = plot.denoised +
  ggtitle("De-noised") +
  theme_void() +
  theme(legend.position="right", legend.title = element_blank(), 
        plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
        axis.title.x=element_text(size=ax_size, family="sans")) +
  xlab(glue(top_gene_name))

ggsave(filename=glue("{directory.images}{name.dataset}_{top_gene_name}_denoised.pdf"), 
       plot=plot.denoised, width=width, height=height, units="in")
print(plot.denoised)


#Repeat for top 10 genes
top_gene_idx <- head(order(scores[ ,1], decreasing = TRUE), 10)
top_gene_name <- rownames(counts.significant)[top_gene_idx]

pt_size = 0.5
ncol=5
height = 15
width = (height*ncol)/2


#Create formatting
manuscript_Scomp_plot_settings <- function(){
  manuscript_plot_settings <- theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(hjust=0.5, size=30, family="sans", face='bold'),
    )
}

denoised_counts.significant <- scores[ ,1:9] %*% t(loadings[ ,1:9])

plot.unnoised <- plot.Scomps(locations=locations.significant, loadings=t(counts.normalized[top_gene_idx, ]),
                              size=pt_size, colormap_continuous=Reds, ncol=ncol,
                              extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix= "gene", is_equal_scale = FALSE)
ggsave(filename=glue("{directory.images}{name.dataset}_top10genes_comp1_unnoised.pdf"), 
        plot=plot.unnoised, width=width, height=height, units="in")
  
plot.denoised <- plot.Scomps(locations=locations.final, loadings=t(denoised_counts.significant[top_gene_idx, ]),
                              size=pt_size, colormap_continuous=Reds, ncol=ncol, 
                              extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix='gene', is_equal_scale=FALSE)
ggsave(filename=glue("{directory.images}{name.dataset}_top10genes_comp1_denoised.pdf"), 
        plot=plot.denoised, width=width, height=height, units="in")






plot.points(locations=locations.significant, t(counts.normalized[top_gene_idx, ]), 
            size=0.1, cmap_continuous=Reds)
plot.points(locations=locations.final, t(denoised_counts.significant[top_gene_idx, ]), 
            size=0.1, cmap_continuous=Reds)


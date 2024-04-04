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

cat.script_title("Plotting HER2 dataset")



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
source("src/src_gene_stats_tests.R")



# |||||||||||||||||||||
# Global variables ----
# |||||||||||||||||||||

cat.section_title("Global variables")

# Directories
directory.initial_data <- "data/HER2/"
directory.results <- "results/HER2/"
directory.images <- "images/HER2/"
directory.images.manuscript <- "images/HER2/manuscript/"
name.dataset <- "HER2"

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



# |||||||||||||||||||||
# Data ----
# |||||||||||||||||||||

# Original Data
load(glue("{directory.initial_data}HER2.RData"))

#SPCA results
load(glue("{directory.results}SPCA_results_HER2.RData"))

#fPCA results
load(glue("{directory.results}HER2_processed.RData"))
load(glue("{directory.results}HER2_mesh.RData"))
load(glue("{directory.results}HER2_mean.RData"))
load(glue("{directory.results}HER2_nComp_selection.RData"))
load(glue("{directory.results}HER2_components.RData"))
load(glue("{directory.results}HER2_clusters_on_HR_grid_nComp3.RData"))



#Defining true labels for fPCA with named atomic vector
true_labels_named_vector <- true_labels$true_label
names(true_labels_named_vector) <- rownames(true_labels)
true_labels_fPCA <- true_labels_named_vector[rownames(locations.final)]



# |||||||||||||||||||||
# Alignment of Clustering Labels ----
# |||||||||||||||||||||

# fPCA LR alignment
fPCA_aligned_list <- get_aligned_labels(clustering_labels = cluster_labels, ground_truth=true_labels_fPCA, return_mapping=TRUE)
fPCA_labels_aligned <- fPCA_aligned_list[["aligned_clustering_labels"]]
fPCA_gt_mapping <- fPCA_aligned_list[["optimal_permutation"]]

#fPCA HR alignment
fPCA.HR_labels_aligned <- align_with_mapping(cluster_labels_HR, fPCA_gt_mapping)


# SPCA alignment: already aligned



# |||||||||||||||||||||
# Scoring of Clustering Labels ----
# |||||||||||||||||||||

#fPCA - LR labels
names.locations <- rownames(locations.final)
names.locations_unknown <- names(true_labels_fPCA[true_labels_fPCA == '7'])
true_labels_no7 <- true_labels_fPCA[!names(true_labels_fPCA) %in% names.locations_unknown]
fPCA_labels_no7 <- fPCA_labels_aligned[!names(fPCA_labels_aligned) %in% names.locations_unknown]

fPCA_ARI <- adjustedRandIndex(fPCA_labels_no7, true_labels_no7)
fPCA_CHAOS <- fx_CHAOS(as.character(true_labels_fPCA), locations.final)

#Calculate scores for HR labels
cluster_labels_HR_subsampled <- c()
for(name.l in names.locations){
  
  distances <- dist_point_from_points(locations.final[name.l,], grid)
  index.closest_point <- which.min(distances)
  cluster_labels_HR_subsampled[name.l] <- cluster_labels_HR[index.closest_point]
  
}

fPCA_labels_HR_subsampled_aligned <- get_aligned_labels(clustering_labels = cluster_labels_HR_subsampled, 
                                                        ground_truth=true_labels_fPCA, return_mapping=FALSE)
fPCA_labels_HR_subsampled_aligned_no7 <- fPCA_labels_HR_subsampled_aligned[!names(fPCA_labels_HR_subsampled_aligned) %in% names.locations_unknown]
  
fPCA_HR_ncomp3_ARI <- adjustedRandIndex(fPCA_labels_HR_subsampled_aligned_no7, true_labels_no7)
fPCA_HR_ncomp3_CHAOS <- fx_CHAOS(fPCA_labels_HR_subsampled_aligned_no7, locations.final)





# |||||||||||||||||||||
# Truth-SPCA-fPCA ----
# |||||||||||||||||||||

cat.subsection_title("Truth-SPCA-fPCA")

if(RUN[["Truth-SPCA-fPCA"]]){
  
  #Plot Ground Truth
  plot_gt <- plot.points(locations=locations, data=true_labels, size=3.5, discrete=TRUE, 
                         cmap_discrete=qual_cmap)
  print(plot_gt)
  
  #Plot Spatial PCA
  plot_SPCA <- plot.points(locations=locations_final_SPCA, data=aligned_SPCA_labels, size=3.5, discrete=TRUE, 
                           cmap_discrete=qual_cmap)
  print(plot_SPCA)
  
  #Plot fPCA
  plot_fPCA <- plot.points(locations=locations.final, data=fPCA_labels_aligned, size=3.5, discrete=TRUE, 
                           cmap_discrete=qual_cmap)
  print(plot_fPCA)
  
  
  if (MAN_RUN){
    
    #Set some plotting variables
    width = 0.86
    height = 0.86
    pt_size = 0.25
    tit_size = 6
    ax_size = 4
    
    #Gound Truth
    plot_gt_man <- plot.points(locations=locations, data=true_labels, size=pt_size, discrete=TRUE, 
                               custom_colors=qual_custom_colors_gt)
    plot_gt_man = plot_gt_man +
      ggtitle("Ground Truth") +
      theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5))
    ggsave(filename=glue("{directory.images}HER2_GroundTruth.pdf"), 
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
    ggsave(filename=glue("{directory.images}HER2_SPCA.pdf"), 
           plot=plot_SPCA_man, width=width, height=height, units="in")
    print(plot_SPCA_man) 
    
    #fPCA
    num_labels_FPCA <- length(unique(fPCA_labels_aligned))
    qual_custom_colors_FPCA <- qual_custom_colors_gt[1:num_labels_FPCA]
    plot_FPCA_man <- plot.points(locations=locations.final, data=fPCA_labels_aligned, size=pt_size, discrete=TRUE, 
                                 custom_colors=qual_custom_colors_FPCA)
    plot_FPCA_man = plot_FPCA_man +
      ggtitle("fPCA") +
      theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
            axis.title.x=element_text(size=ax_size, family="sans")) +
      xlab(glue("ARI: {round(fPCA_ARI,2)}  |  CHAOS: {round(fPCA_CHAOS,2)}"))
    ggsave(filename=glue("{directory.images}HER2_FPCA.pdf"), 
           plot=plot_FPCA_man, width=width, height=height, units="in")
    print(plot_FPCA_man) 
    
  }
  
  
}



# |||||||||||||||||||||
# High-Res fPCA ----
# |||||||||||||||||||||

if(RUN[["High-Res fPCA"]]){
  
  #Set some plotting variables
  width = 2.4
  height = 2.4
  pt_size = 0.2
  tit_size = 16.74
  ax_size = 11.16
  num_labels_FPCA <- length(unique(fPCA_labels_aligned))
  qual_custom_colors_FPCA <- qual_custom_colors_gt[1:num_labels_FPCA]
  
  filepath <- glue("{directory.results}HER2_HRandLR_clustering_scores.RData")
  load(filepath)
  
  plot_FPCA.HR <- plot.points(locations=grid, data=fPCA.HR_labels_aligned, size=0.8, discrete=TRUE, 
                              custom_colors = qual_custom_colors_FPCA)
  print(plot_FPCA.HR)
  
  if (MAN_RUN){
    
    plot_FPCA.HR_man <- plot.points(locations=grid, data=fPCA.HR_labels_aligned, size=pt_size, discrete=TRUE, 
                                    custom_colors = qual_custom_colors_FPCA)
    plot_FPCA.HR_man = plot_FPCA.HR_man +
      ggtitle("fPCA - High Res") +
      theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(size=tit_size, face='bold', family="sans", hjust=0.5),
            axis.title.x=element_text(size=ax_size, family="sans")) +
      xlab(glue("ARI: {round(fPCA_HR_ncomp3_ARI,2)}  |  CHAOS: {round(fPCA_HR_ncomp3_CHAOS,2)}"))
    ggsave(filename=glue("{directory.images}HER2_FPCA_HR.pdf"), 
           plot=plot_FPCA.HR_man, width=width, height=height, units="in")
    print(plot_FPCA.HR_man) 
  }
  
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
            plot.title=element_text(hjust=0.5, size=16, family="sans", face='bold'),
      )
  }
  
  #Plot comps without scaling
  plot_Scomps <- plot.Scomps(locations=locations.final, loadings=loadings, 
                             size=0.25, extra_formatting = manuscript_Scomp_plot_settings,
                             title_prefix="", is_equal_scale=FALSE)
  grid.arrange(plot_Scomps)
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatialComps.pdf"), 
         plot=plot_Scomps, width=4, height=4, units="in")
  
  #Plot comps with scaling
  plot_Scomps_scaled <- plot.Scomps(locations=locations.final, loadings=loadings, 
                                    size=0.25, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=TRUE)
  grid.arrange(plot_Scomps_scaled)
  ggsave(filename=glue("{directory.images}{name.dataset}_FPCA_SpatialComps_Scaled.pdf"), 
         plot=plot_Scomps_scaled, width=4, height=4, units="in")
  
  #Plot comps for manuscript
  if (MAN_RUN){
    #Select comps
    selected_comps <- seq(4)
    FPCA_comp_colors <- coolwarm(200)
    SPCA_comp_colors <- rev(FPCA_comp_colors)
    pt_size = 1.6
    
    #FPCA selected comps without scaling
    plot_Scomps <- plot.Scomps(locations=locations.final, loadings=loadings[ ,selected_comps], 
                               size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                               title_prefix="", is_equal_scale=FALSE, ncol=4, 
                               colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}HER2_FPCA_SpatiaComps1to4_Unscaled.pdf"),
           plot=plot_Scomps, width=8, height=2, units='in')
    grid.arrange(plot_Scomps)
    
    #FPCA selected comps with scaling
    plot_Scomps_scaled <- plot.Scomps(locations=locations.final, loadings=loadings[ ,selected_comps], 
                                      size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                      title_prefix="", is_equal_scale=TRUE, ncol=4,
                                      colormap_continuous = FPCA_comp_colors)
    ggsave(filename=glue("{directory.images}HER2_FPCA_SpatiaComps1to4_Scaled.pdf"),
           plot=plot_Scomps_scaled, width=8, height=2, units='in')
    grid.arrange(plot_Scomps_scaled)
    
    #SPCA selected comps without scaling
    plot_Scomps_SPCA <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_obj@SpatialPCs)[ ,selected_comps],
                               size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                               title_prefix="", is_equal_scale=FALSE, ncol=4,
                               colormap_continuous = SPCA_comp_colors)
    ggsave(filename=glue("{directory.images}HER2_SPCA_SpatiaComps1to4_Unscaled.pdf"),
           plot=plot_Scomps_SPCA, width=8, height=2, units='in')
    grid.arrange(plot_Scomps_SPCA)
    
    #SPCA selected comps with scaling
    plot_Scomps_SPCA_scaled <- plot.Scomps(locations=locations_final_SPCA, loadings=t(SPCA_obj@SpatialPCs)[ ,selected_comps],
                                    size=pt_size, extra_formatting = manuscript_Scomp_plot_settings,
                                    title_prefix="", is_equal_scale=TRUE, ncol=4,
                                    colormap_continuous = SPCA_comp_colors)
    ggsave(filename=glue("{directory.images}HER2_SPCA_SpatiaComps1to4_Scaled.pdf"),
           plot=plot_Scomps_SPCA_scaled, width=8, height=2, units='in')
    grid.arrange(plot_Scomps_SPCA_scaled)
    
  }
  
}




# |||||||||||||||||||||
# ARI/CHAOS Vs. NSCs ----
# |||||||||||||||||||||

cat.subsection_title("ARI/CHAOS Vs. NSCs")

tic()
if(RUN[["ARI/CHAOS Vs. NSCs"]]){
  
  #Load in scoring dataframes
  filepath <- glue("{directory.results}HER2_HRandLR_clustering_scores.RData")
  load(filepath)
  
  #format for plotting
  ARIs_df <- data.frame(LR_scores=score_by_NSCs_LR_df$ARI, HR_scores=score_by_NSCs_HR_df$ARI)
  CHAOS_df <- data.frame(LR_scores=score_by_NSCs_LR_df$CHAOS, HR_scores=score_by_NSCs_HR_df$CHAOS)
  ARIs_long <- rownames_to_column(ARIs_df, var="NSCs")
  CHAOS_long <- rownames_to_column(CHAOS_df, var="NSCs")
  ARIs_long$NSCs <- as.numeric(ARIs_long$NSCs)
  CHAOS_long$NSCs <- as.numeric(CHAOS_long$NSCs)
  ARIs_long <- pivot_longer(ARIs_long, cols= -NSCs, names_to="Clustering Resolution", values_to="Score")
  CHAOS_long <- pivot_longer(CHAOS_long, cols= -NSCs, names_to="Clustering Resolution", values_to="Score")
  ARIs_long$'Clustering Resolution' <- factor(ARIs_long$'Clustering Resolution', 
                                              levels = c("HR_scores", "LR_scores"), 
                                              labels = c("High", "Low"))
  CHAOS_long$'Clustering Resolution' <- factor(CHAOS_long$'Clustering Resolution', 
                                               levels = c("HR_scores", "LR_scores"), 
                                               labels = c("High", "Low"))
  
  #Set plot themes
  plot_themes <- function(){
    plot_themes <- theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.position = c(0.8, 0.15),  # Bottom right
            legend.justification = c(1, 0),  # Anchor point is the bottom right
            legend.box.margin = margin(-5, -5, -5, -5),  # Adjusts the margin to move the legend inside the plot area
            legend.margin = margin(0, 0, 0, 0),
            legend.key.height = unit(0.7,"lines"),
            axis.line = element_line(color = "dimgrey"),
            axis.text.x = element_text(color = "dimgrey", size=8),
            axis.text.y = element_text(color = "dimgrey", size=8),
            axis.title.x = element_text(size=10),
            axis.title.y = element_text(size=10),
            legend.title = element_text(size=9),
            legend.text = element_text(size=7)) 
  }
  
  #Plot ARI vs NSC
  plot.ARIvsNSC <- ggplot(ARIs_long, aes(x=NSCs, y=Score, color=`Clustering Resolution`, group=`Clustering Resolution`)) +
    geom_line() +
    plot_themes() +
    scale_x_continuous(breaks = seq(0, max(ARIs_long$NSCs), by = 2)) +
    scale_color_brewer(palette = qual_cmap) +
    labs(x = "Number of Spatial Components", y = "ARI")
  
  ggsave(filename=glue("{directory.images}HER2_ARIvsNSC.pdf"), 
         plot=plot.ARIvsNSC, width=4, height=2, units="in")
  
  print(plot.ARIvsNSC)
  
  
  #Plot CHAOS vs NSC
  plot.CHAOSvsNSC <- ggplot(CHAOS_long, aes(x=NSCs, y=Score, color=`Clustering Resolution`, group=`Clustering Resolution`)) +
    geom_line() +
    plot_themes() +
    theme(legend.position = c(0.8, 0.7)) +
    scale_x_continuous(breaks = seq(0, max(CHAOS_long$NSCs), by = 2)) +
    scale_color_brewer(palette = qual_cmap) +
    labs(x = "Number of Spatial Components", y = "CHAOS")
  
  ggsave(filename=glue("{directory.images}HER2_CHAOSvsNSC.pdf"), 
         plot=plot.CHAOSvsNSC, width=4, height=2, units="in")
  
  print(plot.CHAOSvsNSC)
  
}


# |||||||||||||||||||||
# Clustering Vs. NSCs ----
# |||||||||||||||||||||


cat.subsection_title("Clustering Vs. NSCs")

if(RUN[["Clustering Vs. NSCs"]]){
  
  #Align all LR clusters with ground truth
  cluster_labels_aligned_LR_all <- as.data.frame(matrix(NA, nrow=dim(cluster_labels_LR_all)[1], 
                                                        ncol=dim(cluster_labels_LR_all)[2]))
  for (i in seq(colnames(cluster_labels_LR_all))){
    LR_labels <- cluster_labels_LR_all[ ,i]
    fPCA_aligned_list <- get_aligned_labels(clustering_labels = LR_labels, 
                                            ground_truth=true_labels_fPCA, return_mapping=TRUE)
    fPCA_labels_aligned <- fPCA_aligned_list[["aligned_clustering_labels"]]
    fPCA_gt_mapping <- fPCA_aligned_list[["optimal_permutation"]]
    cluster_labels_aligned_LR_all[ ,i] <- fPCA_labels_aligned
    
  }
  
  #Downsample HR clusters and align with ground truth
  names.locations <- rownames(locations.final)
  HRtoLR_cluster_labels <- as.data.frame(matrix(NA, nrow=length(names.locations), 
                                                ncol=ncol(cluster_labels_HR_all)))
  rownames(HRtoLR_cluster_labels) <- names.locations
  for(name.l in names.locations){
    
    distances <- dist_point_from_points(locations.final[name.l,], grid)
    index.closest_point <- which.min(distances)
    HRtoLR_cluster_labels[name.l, ] <- cluster_labels_HR_all[index.closest_point, ]
    
  }
  
  cluster_labels_aligned_HR_all <- as.data.frame(matrix(NA, nrow=nrow(HRtoLR_cluster_labels), 
                                                        ncol=ncol(HRtoLR_cluster_labels)))
  for (i in seq(colnames(HRtoLR_cluster_labels))){
    HR_labels <- HRtoLR_cluster_labels[ ,i]
    fPCA_aligned_list <- get_aligned_labels(clustering_labels = HR_labels, 
                                            ground_truth=true_labels_fPCA, return_mapping=TRUE)
    fPCA_labels_aligned <- fPCA_aligned_list[["aligned_clustering_labels"]]
    fPCA_gt_mapping <- fPCA_aligned_list[["optimal_permutation"]]
    cluster_labels_aligned_HR_all[ ,i] <- fPCA_labels_aligned
    
  }
  
  #Create formatting
  manuscript_Scomp_plot_settings <- function(){
    manuscript_plot_settings <- theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(hjust=0.5, size=8, family="sans", face='bold'),
      )
  }
  
  #Define color scheme
  num_labels_FPCA <- length(unique(fPCA_labels_aligned))
  qual_custom_colors_FPCA <- qual_custom_colors_gt[1:num_labels_FPCA]
  
  #Plot all clusters - LR
  plot.clusterings_by_NSCs_all_LR <- plot.Scomps(locations.final, loadings=cluster_labels_aligned_LR_all, 
                                                 size=0.25, extra_formatting = manuscript_Scomp_plot_settings,
                                                 title_prefix="", is_equal_scale=FALSE, colormap_discrete=qual_cmap, 
                                                 discrete=TRUE, custom_colors = qual_custom_colors_FPCA)
  ggsave(filename=glue("{directory.images}HER2_AllClusterings_By_NSCs_LR.pdf"), 
         plot=plot.clusterings_by_NSCs_all_LR, width=4, height=4, units="in")
  grid.draw(plot.clusterings_by_NSCs_all_LR)
  
  
  #Plot all clusters - HR
  plot.clusterings_by_NSCs_all_HR <- plot.Scomps(locations.final, loadings=cluster_labels_aligned_HR_all, 
                                                 size=0.25, extra_formatting = manuscript_Scomp_plot_settings,
                                                 title_prefix="", is_equal_scale=FALSE, colormap_discrete=qual_cmap,
                                                 discrete=TRUE, custom_colors=qual_custom_colors_FPCA)
  ggsave(filename=glue("{directory.images}HER2_AllClusterings_By_NSCs_HRtoLR.pdf"), 
         plot=plot.clusterings_by_NSCs_all_HR, width=4, height=4, units="in")
  grid.draw(plot.clusterings_by_NSCs_all_HR)
  
  
  #Select specific clusters - LR
  selected_idx_clusterings <- c(1, 2, 3, 4)
  selected_LR_clusterings <- as.data.frame(matrix(NA, nrow=dim(cluster_labels_LR_all)[1], 
                                                  ncol=length(selected_idx_clusterings)))
 
  for (i in seq(length(selected_idx_clusterings))){
    idx <- selected_idx_clusterings[i]
    LR_clusters <- cluster_labels_aligned_LR_all[ ,idx]
    selected_LR_clusterings[ ,i] <- LR_clusters
  }
  
  #Select specific clusters - HR to LR
  selected_idx_clusterings <- c(1, 2, 3, 4)
  selected_HRtoLR_clusterings <- as.data.frame(matrix(NA, nrow=dim(cluster_labels_aligned_HR_all)[1], 
                                                  ncol=length(selected_idx_clusterings)))
  
  for (i in seq(length(selected_idx_clusterings))){
    idx <- selected_idx_clusterings[i]
    HRtoLR_clusters <- cluster_labels_aligned_HR_all[ ,idx]
    selected_HRtoLR_clusterings[ ,i] <- HRtoLR_clusters
  }
  
  #Create formatting
  manuscript_selected_clusters_settings <- function(){
    manuscript_plot_settings <- theme_void() +
      theme(legend.position="none", 
            plot.title=element_text(hjust=0.5, size=16, family="sans", face='bold'),
      )
  }

  #Plot specific clusters - LR
  custom_titles <- c("1", "1-2", "1-3", "1-4")
  plot.selected_clusters_LR_by_NSC <- plot.Scomps(locations.final, loadings=selected_LR_clusterings, 
                                                  size=1.5, extra_formatting = manuscript_selected_clusters_settings,
                                                  title_prefix="", is_equal_scale=FALSE, colormap_discrete=qual_cmap, 
                                                  discrete=TRUE, ncol=4, custom_colors = qual_custom_colors_FPCA,
                                                  custom_titles = custom_titles)
  ggsave(filename=glue("{directory.images}HER2_SelectedClusterings_By_NSCs_LR.pdf"), 
         plot=plot.selected_clusters_LR_by_NSC, width=8, height=2, units="in")
  grid.draw(plot.selected_clusters_LR_by_NSC)
  
  
  #Plot specific clusters - HR to LR
  plot.selected_clusters_HRtoLR_by_NSC <- plot.Scomps(locations.final, loadings=selected_HRtoLR_clusterings, 
                                                  size=1.5, extra_formatting = manuscript_selected_clusters_settings,
                                                  title_prefix="", is_equal_scale=FALSE, colormap_discrete=qual_cmap, 
                                                  discrete=TRUE, ncol=4, custom_colors = qual_custom_colors_FPCA,
                                                  custom_titles = custom_titles)
  ggsave(filename=glue("{directory.images}HER2_SelectedClusterings_By_NSCs_HRtoLR.pdf"), 
         plot=plot.selected_clusters_HRtoLR_by_NSC, width=8, height=2, units="in")
  grid.draw(plot.selected_clusters_HRtoLR_by_NSC)
  
  
}




# |||||||||||||||||||||
# Gene by Cluster ----
# |||||||||||||||||||||

cat.subsection_title("Gene by Cluster")

tic()
if(RUN[["Gene by Cluster"]]){
  
  #Select cancer labels
  labels_of_interest <- c(3,6)
  
  #Run 2 sample t-test
  p_adj_values <- between_cluster_ttest(fPCA_labels_aligned, labels_of_interest, counts.significant)
  
  #Get significant genes
  sig_genes <- names(p_adj_values[p_adj_values < 1e-10])
  num_sig_genes <- length(sig_genes)
  top_n_genes <- head(sort(sig_genes,decreasing=TRUE), 25)
  
  #Get top n fold_change sig genes
  top_genes_by_cluster_by_df <- get_gene_means_by_cluster(fPCA_labels_aligned, counts.significant, sig_genes)
  top_genes_by_cluster_by_df$fold_change <- apply(top_genes_by_cluster_by_df[ ,labels_of_interest], 1, mean) / apply(top_genes_by_cluster_by_df[ ,-labels_of_interest], 1, mean)
  top_genes_by_cluster_by_df <- top_genes_by_cluster_by_df %>% arrange(desc(fold_change))
  top_genes_by_cluster_by_df <- top_genes_by_cluster_by_df[1:9, 1:6]
  
  #Z-score all genes for easier visualization and comparison
  z_score_df <- apply(top_genes_by_cluster_by_df, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  z_score_df <- t(z_score_df)
  z_score_df <- as.data.frame(z_score_df)
  
  #Plot significant genes
  top_genes_by_cluster_by_df_melt  <- tibble::rownames_to_column(z_score_df, "Gene Name")
  top_genes_by_cluster_by_df_melt  <- melt(top_genes_by_cluster_by_df_melt, id.vars="Gene Name", variable.name="Cluster", value.name="Gene Mean")
  
  ggplot(top_genes_by_cluster_by_df_melt, aes(x=Cluster, y=`Gene Mean`, fill=factor(`Gene Name`))) +
    geom_bar(stat='identity', position='dodge') +
    theme(legend.position="bottom") +
    scale_fill_brewer(palette="Spectral") +
    ylab("Gene Mean Z-scores")
  
  
  if (MAN_RUN){
    
    #Crete custom selected gene df for plotting
    selected_genes <- c("ERBB2", "MGP", "MDK", "PEG10", "SOX11", "CFD")
    sel_genes_by_cluster_df <- get_gene_means_by_cluster(fPCA_labels_aligned, counts.significant, selected_genes)
    z_score_df <- apply(sel_genes_by_cluster_df, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    z_score_df <- as.data.frame(t(z_score_df))
    genecluster_df <- tibble::rownames_to_column(z_score_df, "Gene Name")
    col_labels <- c("Gene Name", "Fat Tissue", "Normal Glands", "DCIS", "Fibrous Tissue", "Immune Cells","Invasive Cancer")
    colnames(genecluster_df) <- col_labels
    genecluster_df_melt  <- melt(genecluster_df, id.vars="Gene Name", variable.name="Cluster", value.name="Gene Mean")
    
    #Reorganize so I can plot cluster labels in specific order
    cluster_order <- c("DCIS", "Invasive Cancer", "Immune Cells", "Fibrous Tissue","Normal Glands", "Fat Tissue")
    genecluster_df_melt$Cluster <- factor(genecluster_df_melt$Cluster, levels=cluster_order)
    
    #Plot 
    plot.gene_by_cluster <- ggplot(genecluster_df_melt, aes(x=Cluster, y=`Gene Mean`, fill=factor(`Gene Name`))) +
      geom_bar(stat='identity', position='dodge') +
      theme(legend.position='right', 
            axis.text.x = element_text(angle=45, hjust=1, color = "dimgrey", size=10),
            axis.text.y = element_text(color = "dimgrey", size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=14, 
                                        margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +
      scale_fill_brewer(palette="Set2") +
      ylab("Gene Mean Z-Scores") +
      labs(fill="Gene")
    
    ggsave(filename=glue("{directory.images}GeneMean_by_cluster_cancer.pdf"), 
           plot=plot.gene_by_cluster, width=6, height=3.2, units="in")
    print(plot.gene_by_cluster)
    
    
  }
  
}



# |||||||||||||||||||||
# Volcano Plot ----
# |||||||||||||||||||||

# Select cancer labels and organize data by q value and fold change
labels_of_interest <- c(3,6) #Cancer Labels
genes_by_cluster_df <- get_gene_means_by_cluster(fPCA_labels_aligned, counts.significant, 
                                                 names(p_adj_values))
genes_by_cluster_df$log2_fold_change <- log2(apply(genes_by_cluster_df[ ,labels_of_interest], 1, mean) / 
                                               apply(genes_by_cluster_df[ ,-labels_of_interest], 1, mean))
volcano_df <- data.frame(gene_name=rownames(genes_by_cluster_df), q_val=p_adj_values, log2_fold_change=genes_by_cluster_df$log2_fold_change)

# Select specific labels for plotting
volc_df_filt <- volcano_df[volcano_df$q_val < 1e-10, ]
top_n_sig <- rownames(head(volc_df_filt[order(volc_df_filt$q_val), ], 10))
top_n_FC <- rownames(head(volc_df_filt[order(abs(volc_df_filt$log2_fold_change), decreasing=TRUE), ], 9))
lab_select <- union(top_n_sig, top_n_FC)

# Plot volcano plot
plot.volcano_cancer <- EnhancedVolcano(volcano_df, lab=volcano_df$gene_name, x="log2_fold_change", y="q_val", drawConnectors = TRUE,
                widthConnectors = 0.5, pCutoff = 1e-10, FCcutoff = 1.0, max.overlaps = Inf, labSize=2.5,
                caption = bquote(""), title="", subtitle="", selectLab = lab_select) +
                #legendPosition = "right", legendLabSize = 8)
  theme(legend.position='none', axis.line = element_line(color = "dimgrey"),
        axis.text.x = element_text(color = "dimgrey", size=10),
        axis.text.y = element_text(color = "dimgrey", size=10),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  coord_flip()

ggsave(filename=glue("{directory.images}HER2_VolcanoPlot_CancerVsNoncancer.pdf"), 
       plot=plot.volcano_cancer, width=6, height=4, units="in")
print(plot.volcano_cancer)


###Plot cluster comparisons
#Set some plotting variables
scale_factor=1
width = 0.86*scale_factor
height = 0.86*scale_factor
pt_size = 0.25
tit_size = 6*scale_factor
ax_size = 4*scale_factor

#Cancer Cluster
cancer_mapping <- t(rbind(c(1,2,3,4,5,6,7), c(1,1,2,1,1,2,1)))
cancer_colors <- c("#BEAED4", "#F0027F")
cancer_labels <- align_with_nonunique_mapping(clustering_labels=fPCA_labels_aligned, 
                                    mapping = cancer_mapping)
plot_cancer_cluster <- plot.points(locations=locations.final, data=cancer_labels, size=pt_size, discrete=TRUE, 
                           custom_colors=cancer_colors)
plot_cancer_cluster = plot_cancer_cluster +
  ggtitle("Ground Truth") +
  theme_void() +
  theme(legend.position="none", 
        plot.title=element_blank())
ggsave(filename=glue("{directory.images}HER2_CancerCluster_highlighted.pdf"), 
       plot=plot_cancer_cluster, width=width, height=height, units="in")
print(plot_cancer_cluster) 

  

# |||||||||||||||||||||
# Hierarchical Clustering by Cluster Labels ----
# |||||||||||||||||||||
#Organize Gene mean by cluster label matrix
gene_means_by_cluster <- get_gene_means_by_cluster(fPCA_labels_aligned, counts.significant, rownames(counts.significant))
gene_zscores_by_cluster <- t( apply(gene_means_by_cluster, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)) )                

#Get label names for cluster number and genes
gt_labels <- c("Fat Tissue", "Normal Glands", "DCIS", "Fibrous Tissue", "Immune Cells","Invasive Cancer")
colnames(gene_zscores_by_cluster) <- gt_labels
selected_genes <- c("ERBB2", "MDK", "PEG10", "SOX11", "CFD")
no_gene_labels <- as.list(rep("", length(rownames(gene_zscores_by_cluster))))
selected_genes_full <- as.list(rep("", length(rownames(gene_zscores_by_cluster))))
for (gene in selected_genes){
  idx <- which(rownames(gene_zscores_by_cluster) == gene)
  selected_genes_full[idx] <- gene
}

#Make heatmap
scale_factor = 3
width = 2.32*scale_factor
height = 1.3*scale_factor
plot.heatmap <- Heatmap(matrix=t(gene_zscores_by_cluster), name='z-score',
                        row_names_side='left', row_dend_side='right', column_names_side='top',
                        column_dend_side='bottom', column_labels = no_gene_labels,
                        column_names_gp=gpar(fontsize=5), row_names_gp=gpar(fontsize=15), 
                        column_title='Genes', column_title_gp = gpar(fontsize=15))
print(plot.heatmap)
pdf(glue("{directory.images}{name.dataset}_genebycluster_heatmap_test.pdf"), width=width, height=height)
plot.heatmap
dev.off()










###Pietro's test code for clusteirng H matrix
scores_new <- scores
for(h in 1:3){
  scores_new[,h] <- scores_new[,h]*sqrt(sum(loadings[,h]^2))
}
plot3d(scores_new[,1:3], col = color)


# Select cancer labels and organize data by q value and fold change
labels_of_interest <- c(3,6) #Cancer Labels
genes_by_cluster_df <- get_gene_means_by_cluster(fPCA_labels_aligned, counts.significant, 
                                                 names(p_adj_values))
genes_by_cluster_df$log2_fold_change <- log2(apply(genes_by_cluster_df[ ,labels_of_interest], 1, mean) / 
                                               apply(genes_by_cluster_df[ ,-labels_of_interest], 1, mean))
volcano_df <- data.frame(gene_name=rownames(genes_by_cluster_df), q_val=p_adj_values, log2_fold_change=genes_by_cluster_df$log2_fold_change)

# Select specific labels for plotting
volc_df_filt <- volcano_df[volcano_df$q_val < 1e-10 &volcano_df$log2_fold_change > 1, ]



indexes <- which(volcano_df$q_val < 1e-10 &volcano_df$log2_fold_change > 1)
color <- rep("black", nrow(scores))
color[indexes] <- "red"

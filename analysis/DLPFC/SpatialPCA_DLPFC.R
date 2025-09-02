## ||||||||||||||||||||
# Install and load in libraries ----
## ||||||||||||||||||||

### Installation of packages from github
library(devtools)
if (!requireNamespace("SPARK", quietly=TRUE)) {devtools::install_github('xzhoulab/SPARK')}
if (!requireNamespace("SpatialPCA", quietly=TRUE)) {devtools::install_github("shangll123/SpatialPCA")}

### Load in all libraries
library(glue)
library(SPARK)
library(Seurat)
library(peakRAM)
library(SpatialPCA)
library(ggplot2)
library(tictoc)
library(dplyr)
tic.clearlog()





#=============================================================================#




## ||||||||||||||||||||
# Define filepaths for loading in and saving data ----
## ||||||||||||||||||||

###Directories for loading/saving in files/data (REPLACE WITH PERSONAL PATHS)
project_dir <- "/Users/aburns/Codebook/Projects/spatial_transcriptomics_GRSVD/Paper/ST-lens_NEW/"
source_dir <- glue("{project_dir}src/")
data_dir <- glue("{project_dir}data/DLPFC/")
results_dir <- glue("{project_dir}results/DLPFC/")
DLPFC_sample_9_filepath <- glue("{data_dir}data.RData")
name.dataset <- "DLPFC_sample9"

### Load in functions/variables form source (src) files
source(glue("{source_dir}src_support_fxns.R"))

### Define run logic gates
RUN_time_trials <- F




#=============================================================================#




## ||||||||||||||||||||
# Load in data for DLPFC sample 9 (151673) ----
## ||||||||||||||||||||

### Load in sample 9
load(DLPFC_sample_9_filepath) #Loads in count_sub, KRM_manual_layers_sub, and xy_coords
xy_coords <- as.matrix(locations)
rownames(xy_coords) <- colnames(counts) #Necessary for SpatialPCA functions




#=============================================================================#




## ||||||||||||||||||||
# Run SpatialPCA (Only need to be run once. Skip to next section if already run) ----
## ||||||||||||||||||||

### Create SpatialPCA object
DLPFC_SPCA <- CreateSpatialPCAObject(counts=counts, location=xy_coords, 
                                     project='SpatialPCA', gene.type="spatial", 
                                     sparkversion="sparkx", numCores_spark=5, 
                                     gene.number=3000, customGenelist=NULL, 
                                     min.loctions = 20, min.features=20)

### Run rest of SpatialPCA pipeline
mem <- peakRAM({
  tic("Overall time")
  tic("Building Kernel")
  DLPFC_SPCA <- SpatialPCA_buildKernel(DLPFC_SPCA, kerneltype='gaussian', bandwidthtype="SJ", bandwidth.set.by.user=NULL)
  toc() #Estimate time: 2 sec (mac M1 chip)
  tic("Estimate Loading")
  DLPFC_SPCA <- SpatialPCA_EstimateLoading(DLPFC_SPCA, fast=FALSE, SpatialPCnum=20)
  toc() #Estimate time: 27.5 min - 45 min (mac M1 chip)
  tic("Calculating Spatial PCs")
  DLPFC_SPCA <- SpatialPCA_SpatialPCs(DLPFC_SPCA, fast=FALSE)
  toc() #Estimate time: 1.3 min (mac M1 chip)
  toc() #Estimate time: 29 min (mac M1 chip)
})




#=============================================================================#




## ||||||||||||||||||||
# TIME TRIALS: Run Spatial PCA ----
## ||||||||||||||||||||

if (RUN_time_trials){
  n_trials <- 4
  trial_list <- list()
  
  for (i in 1:n_trials){
    ### Clear timing log and print iteration number
    tic.clearlog()
    cat("Iteration:", i, "\n")
    
    tic("Overall time")
    ### Create SpatialPCA object
    tic("Create SPCA object")
    SPCA_obj_timetrial <- CreateSpatialPCAObject(counts=counts, location=xy_coords, 
                                         project='SpatialPCA', gene.type="spatial", 
                                         sparkversion="sparkx", numCores_spark=5, 
                                         gene.number=3000, customGenelist=NULL, 
                                         min.loctions = 20, min.features=20)
    toc(log=T) 
    
    ###Run rest of SpatialPCA pipeline
    mem <- peakRAM({
      tic("Building Kernel") 
      SPCA_obj_timetrial = SpatialPCA_buildKernel(SPCA_obj_timetrial, kerneltype="gaussian", bandwidthtype="SJ")
      toc(log=T) 
      tic("Estimate Loading")
      SPCA_obj_timetrial = SpatialPCA_EstimateLoading(SPCA_obj_timetrial,fast=FALSE,SpatialPCnum=20)
      toc(log=T)
      tic("Calculating Spatial PCs")
      SPCA_obj_timetrial = SpatialPCA_SpatialPCs(SPCA_obj_timetrial, fast=FALSE)
      toc(log=T) 
    })
    toc(log=T) #Estimate time: 41 sec (mac M1 chip)
    
    ### Store time trials
    times_log <- tic.log(format=F)
    trial_df <- data.frame(
      Iteration = i,
      Step = sapply(times_log, function(t) t$msg),
      ElapsedTime = sapply(times_log, function(t) t$toc[[1]] - t$tic[[1]])
    )
    trial_list[[i]] <- trial_df
    
    rm(SPCA_obj_timetrial)
    
  }
  
  SPCA_timetrials_df <- bind_rows(trial_list)
  
  ##Save time trials
  filepath <- glue("{results_dir}SPCA_timetrialresults_{name.dataset}.RData")
  save(SPCA_timetrials_df,
       file=filepath)
  
}




#=============================================================================#




## ||||||||||||||||||||
# Walktrap CLustering Algorithm ----
## ||||||||||||||||||||

### Run clustering
clusternum = 7 #Based on SpatialPCA's reported parameters
knn = 70 #Based on SpatialPCA's reported parameters
SPCA_clusterlabels <- walktrap_clustering(clusternum = clusternum, latent_dat=DLPFC_SPCA@SpatialPCs, knearest = knn)
SPCA_clusterlabels_refined <- refine_cluster_10x(clusterlabels = SPCA_clusterlabels, location=DLPFC_SPCA@location, shape="hexagon")




#=============================================================================#




## ||||||||||||||||||||
# Extract and save ground truth labels ----
## ||||||||||||||||||||

### Extract true labels as strings
#old_colnames <- colnames(DLPFC_SPCA@normalized_expr)
#new_colnames <- gsub(".1", "-1", old_colnames)
#colnames(DLPFC_SPCA@normalized_expr) <- new_colnames
#true_labels_str <- KRM_manual_layers_sub$layer_guess_reordered[match(colnames(DLPFC_SPCA@normalized_expr),colnames(count_sub))]

### Convert true labels to integers
true_labels_str <- as.character(true_labels$true_label)
true_labels_int_vec <- as.integer(true_labels$true_label)

### Save true labels with NAs
true_labels_df <- data.frame('anat_labels'=true_labels_str, 'int_labels'=true_labels_int_vec)
rownames(true_labels_df) <- rownames(xy_coords)
true_labels_noNA_df <- true_labels_df[!is.na(true_labels_df[,2]), ]




#=============================================================================#




## ||||||||||||||||||||
# Redo labels based off ARI mapping (Mapped to ground truth) ----
## ||||||||||||||||||||

###Remove label NAs from both true labels and Spatial PCA labels and xy_coords
label_NAs <- which(is.na(true_labels_df[,1]))
SPCA_clusterlabels_refined_noNAs <- SPCA_clusterlabels_refined[-label_NAs]
xy_coords_noNAs <- xy_coords[-label_NAs, ]

###Convert true labels into characters
true_labels_noNA_df$intstr_labels <- as.character(true_labels_noNA_df$int_labels)

###Aligned labels
aligned_SPCA_labels_noNA <- get_aligned_labels(clustering_labels=SPCA_clusterlabels_refined_noNAs, ground_truth=true_labels_noNA_df$intstr_labels)




#=============================================================================#




## ||||||||||||||||||||
# Calculate Adjusted Rand Index and CHAOS ----
## ||||||||||||||||||||

SPCA_ARI <- adjustedRandIndex(aligned_SPCA_labels_noNA, true_labels_noNA_df$intstr_labels)
SPCA_CHAOS <- fx_CHAOS(aligned_SPCA_labels_noNA, xy_coords_noNAs)




#=============================================================================#




## ||||||||||||||||||||
# Plot clustering labels of SPCA and Ground Truth ----
## ||||||||||||||||||||

###Define colors
cbp=c("#5CB85C" ,"#9C9EDE" ,"#FFDC91", "#4DBBD5" ,"#FF9896" ,"#FED439", "#E377C2", "#FED439")

###Ground Truth
plot_cluster(location=xy_coords_noNAs, true_labels_noNA_df$intstr_labels, pointsize=1.5, title_in=paste0("Ground truth"), color_in=cbp)

###SPCA
plot_cluster(location=xy_coords_noNAs, aligned_SPCA_labels_noNA, pointsize=1.5, title_in=paste0("Spatial-PCA"), color_in=cbp)




#=============================================================================#




## ||||||||||||||||||||
## Save (or load) files ----
## ||||||||||||||||||||
###Save all objects in one Rdata file
name.dataset <- "DLPFC_sample9"
SPCA_obj <- DLPFC_SPCA
locations_final_SPCA <- xy_coords_noNAs
aligned_SPCA_labels <- aligned_SPCA_labels_noNA
filepath <- glue("{results_dir}SPCA_allresults_NEW_{name.dataset}.RData")
save(SPCA_obj, xy_coords, locations_final_SPCA, aligned_SPCA_labels, SPCA_ARI, SPCA_CHAOS, file=filepath)








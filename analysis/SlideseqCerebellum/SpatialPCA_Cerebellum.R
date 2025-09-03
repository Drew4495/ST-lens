## ||||||||||||||||||||
# Install and load in libraries ----
## ||||||||||||||||||||

###Installation of packages from github
library(devtools)
if (!requireNamespace("SPARK", quietly=TRUE)) {devtools::install_github('xzhoulab/SPARK')}
if (!requireNamespace("SpatialPCA", quietly=TRUE)) {devtools::install_github("shangll123/SpatialPCA")}

###Load in all libraries
library(SPARK)
library(Seurat)
library(peakRAM)
library(SpatialPCA)
library(ggplot2)
library(glue)
library(tictoc)
library(dplyr)
tic.clearlog()




#=============================================================================#




## ||||||||||||||||||||
# Define filepaths for loading in and saving data ----
## ||||||||||||||||||||

###Directories for loading/saving in files/data (REPLACE WITH PERSONAL PATHS)
project_dir <- paste0(getwd(), "/")
source_dir <- glue("{project_dir}src/")
data_dir <- glue("{project_dir}data/SlideseqCerebellum/")
results_dir <- glue("{project_dir}results/SlideseqCerebellum/")
cerebellum_data_filepath <- glue("{data_dir}Cerebellum_data_for_SPCA/mouse_cerebellum_slideseq_data_lulushang.rds")
pseudo_ground_truth_labels_filepath <- glue("{data_dir}Cerebellum_data_for_SPCA/Cerebellum_PseudoTruth_lbls_lambd1eneg3_knn220_npc20.csv")
pseudotruth_coords_filepath <- glue("{data_dir}Cerebellum_data_for_SPCA/Cerebellum_PostFPCA_locations_PIETRO.csv")
name.dataset <- "SlideseqCerebellum"

###Load in functions/variables form source (src) files
source(glue("{source_dir}src_support_fxns.R"))

### Define run logic gates
RUN_time_trials <- F




#=============================================================================#




## ||||||||||||||||||||
# Load in data ----
## ||||||||||||||||||||

load(cerebellum_data_filepath) #loads in two objects [location and sp_count]
print(dim(sp_count)) # The count matrix (17,729genes x 25,551spots) 
print(dim(location)) # The location matrix (25,551spots x 2)

###Load pseudotruth labels and coords for alignment
pseudotruth_lbls <- read.csv(pseudo_ground_truth_labels_filepath, header=TRUE)
pseudotruth_coords <- read.csv(pseudotruth_coords_filepath, header=TRUE)




#=============================================================================#




## ||||||||||||||||||||
# Run Spatial PCA ----
## ||||||||||||||||||||

tic("Overall time")
###Create SpatialPCA object
tic("Create SPCA object")
SPCA_obj = CreateSpatialPCAObject(counts=sp_count, location=location, project = "SpatialPCA",
                                  gene.type="spatial", sparkversion="sparkx", numCores_spark=5, 
                                  customGenelist=NULL,min.loctions = 20, min.features=20)
toc()

###Run rest of SpatialPCA pipeline
mem <- peakRAM({
  tic("Building Kernel")
  SPCA_obj = SpatialPCA_buildKernel(SPCA_obj, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
  toc() #Estimate time: 69 sec (mac M1 chip) 
  tic("Estimate Loading")
  SPCA_obj = SpatialPCA_EstimateLoading(SPCA_obj,fast=TRUE,SpatialPCnum=20)
  toc() #Estimate time: 1300 sec (mac M1 chip)
  tic("Calculating Spatial PCs")
  SPCA_obj = SpatialPCA_SpatialPCs(SPCA_obj, fast=TRUE)
  toc() #Estimate time: 764 sec (mac M1 chip)
})
toc() #Estimate time: 2400 sec (mac M1 chip)

###Save new locations
locations_postSPCA <- SPCA_obj@location




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
    SPCA_obj_timetrial = CreateSpatialPCAObject(counts=sp_count, location=location, project = "SpatialPCA",
                                      gene.type="spatial", sparkversion="sparkx", numCores_spark=5, 
                                      customGenelist=NULL,min.loctions = 20, min.features=20)
    toc(log=T) 
    
    ###Run rest of SpatialPCA pipeline
    mem <- peakRAM({
      tic("Building Kernel") 
      SPCA_obj_timetrial = SpatialPCA_buildKernel(SPCA_obj_timetrial, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
      toc(log=T) 
      tic("Estimate Loading")
      SPCA_obj_timetrial = SpatialPCA_EstimateLoading(SPCA_obj_timetrial,fast=TRUE,SpatialPCnum=20)
      toc(log=T)
      tic("Calculating Spatial PCs")
      SPCA_obj_timetrial = SpatialPCA_SpatialPCs(SPCA_obj_timetrial, fast=TRUE)
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
# Louvain CLustering Algorithm ----
## ||||||||||||||||||||

###Run clustering
tic("clustering")
clusternum <- 8
knn <- round(sqrt(dim(SPCA_obj@SpatialPCs)[2]))
SpatialPCs <- as.matrix(SPCA_obj@SpatialPCs)
SPCA_8clusterlabels <- louvain_clustering(clusternum = clusternum, latent_dat=SpatialPCs, knearest = knn)
toc()

###With 9 clustering
tic("clustering")
clusternum <- 9
knn <- round(sqrt(dim(SPCA_obj@SpatialPCs)[2]))
SpatialPCs <- as.matrix(SPCA_obj@SpatialPCs)
SPCA_9clusterlabels <- louvain_clustering(clusternum = clusternum, latent_dat=SpatialPCs, knearest = knn)
toc()




#=============================================================================#




## ||||||||||||||||||||
# Align labels based on pseudo ground truth ----
## ||||||||||||||||||||

### Reformat coords for SPCA and pseudotruth
colnames(pseudotruth_coords) <- colnames(locations_postSPCA)

###Find indices in pseudotruth labels that aren't ins SPCA and remove from pseudotruth
indices_no_match_remove <- which(!(pseudotruth_coords$xcoord %in% locations_postSPCA[,1] & pseudotruth_coords$ycoord %in% locations_postSPCA[,2]))
pseudotruth_lbls_matched <- as.matrix(pseudotruth_lbls)[-indices_no_match_remove]
pseudotruth_coords_filtered <- pseudotruth_coords[-indices_no_match_remove,]

###Find indices in SPCA  that's arent in pseudotruth 
indices_no_match_add <- which(!(locations_postSPCA[,1] %in% pseudotruth_coords_filtered$xcoord & locations_postSPCA[,2] %in% pseudotruth_coords_filtered$ycoord))

###Fill in pseudo_truth_lbls (FPCA) with random labels (1-8) that are missing compared to SPCA_coords
#pseudotruth_lbls_matched <- as.matrix(pseudotruth_lbls)
for (idx in indices_no_match_add){
  rand_val <- as.character(sample(1:8, 1, replace = TRUE))
  pseudotruth_lbls_matched <- c(pseudotruth_lbls_matched[1:(idx-1)], rand_val, 
                                pseudotruth_lbls_matched[idx:length(pseudotruth_lbls_matched)])
}

###Redo labels
aligned_SPCA_labels_cluster8 <- get_aligned_labels(clustering_labels=SPCA_8clusterlabels, ground_truth=pseudotruth_lbls_matched)
aligned_SPCA_labels_cluster9 <- get_aligned_labels(clustering_labels=SPCA_9clusterlabels, ground_truth=pseudotruth_lbls_matched)





#=============================================================================#




## ||||||||||||||||||||
# Calculate PAS and CHAOS scores and save ----
## ||||||||||||||||||||

###Calculate CHAOS and PAS
CHAOS_score_8 <- fx_CHAOS(aligned_SPCA_labels_cluster8, locations_postSPCA)
PAS_score_8 <- fx_PAS(aligned_SPCA_labels_cluster8, locations_postSPCA)
CHAOS_PAS_df_8 <- data.frame(CHAOS=CHAOS_score_8, PAS=PAS_score_8)

CHAOS_score_9 <- fx_CHAOS(aligned_SPCA_labels_cluster9, locations_postSPCA)
PAS_score_9 <- fx_PAS(aligned_SPCA_labels_cluster9, locations_postSPCA)
CHAOS_PAS_df_9 <- data.frame(CHAOS=CHAOS_score_9, PAS=PAS_score_9)




## ||||||||||||||||||||
# Plot spatial components ----
## ||||||||||||||||||||
###Define colors
cbp=c("#F5F5DC","#FF7070" , "#C3FFFF" ,"#00E1A1","#256BAA","#FED439" ,"#4DBBD5", "#87CEEB")

###Pseudo-Ground_truth
plot_cluster(location=pseudotruth_coords, as.character(pseudotruth_lbls[,1]), pointsize=1.5,
             title_in="Pseudotruth Labels",color_in=cbp,legend="right")

###SpatialPCA
plot_cluster(location=locations_postSPCA, SPCA_8clusterlabels, pointsize=1.5,
             title_in=glue("SpatialPCA (CHAOS: {round(CHAOS_score_8,3)}, PAS: {round(PAS_score_8,3)})"), color_in=cbp,legend="right")

plot_cluster(location=locations_postSPCA, SPCA_9clusterlabels, pointsize=1.5,
             title_in=glue("SpatialPCA (CHAOS: {round(CHAOS_score_9,3)}, PAS: {round(PAS_score_9,3)})"), color_in=append(cbp,"#A9A9A9",after=4),legend="right")

plot_cluster(location=locations_postSPCA, aligned_SPCA_labels_cluster8, pointsize=1.5,
             title_in=glue("SpatialPCA (CHAOS: {round(CHAOS_score_8,3)}, PAS: {round(PAS_score_8,3)})"), color_in=cbp,legend="right")

plot_cluster(location=locations_postSPCA, aligned_SPCA_labels_cluster9, pointsize=1.5,
             title_in=glue("SpatialPCA (CHAOS: {round(CHAOS_score_9,3)}, PAS: {round(PAS_score_9,3)})"), color_in=append(cbp,"#A9A9A9"),legend="right")




#=============================================================================#




## ||||||||||||||||||||
# Save/Load files ----
## ||||||||||||||||||||

###Save all objects in one Rdata file
locations_preSPCA <- location
locations_final_SPCA <- locations_postSPCA
SPCA_CHAOS_8 <- CHAOS_score_8
SPCA_PAS_8 <- PAS_score_8
SPCA_CHAOS_9 <- CHAOS_score_9
SPCA_PAS_9 <- PAS_score_9
filepath <- glue("{results_dir}SPCA_results_{name.dataset}.RData")
save(SPCA_obj, locations_preSPCA, locations_final_SPCA, 
     SPCA_PAS_8, SPCA_CHAOS_8, SPCA_PAS_9, SPCA_CHAOS_9,
     aligned_SPCA_labels_cluster8, aligned_SPCA_labels_cluster9, 
     SPCA_8clusterlabels, SPCA_9clusterlabels, file=filepath)





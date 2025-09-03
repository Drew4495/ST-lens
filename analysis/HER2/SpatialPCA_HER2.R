## ||||||||||||||||||||
## Install and load in libraries ----
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

### Define run logic gates
RUN_time_trials <- F




## ||||||||||||||||||||
## Define filepaths for loading in and saving data ----
## ||||||||||||||||||||

###Directories for loading/saving in files/data (REPLACE WITH PERSONAL PATHS)
setwd(getwd())
project_dir <- paste0(getwd(), "/") 
source_dir <- glue("{project_dir}src/")
data_dir <- glue("{project_dir}data/HER2/")
results_dir <- glue("{project_dir}results/HER2/")
HER2_SPCA_preprocessed_filepath <- glue("{data_dir}HER2_data_for_SpatialPCA/Tumor_data.rdata")
HER2_H1_raw_filepath <- glue("{data_dir}HER2_data_for_SpatialPCA/H1_selection.tsv")
HER2_gt_labels <- glue("{data_dir}HER2_data_for_SpatialPCA/H1_labeled_coordinates_GROUNDTRUTH.tsv")

###Load in functions/variables form source (src) files
source(glue("{source_dir}src_support_fxns.R"))




## ||||||||||||||||||||
## Load in data for HER2 sample ----
## ||||||||||||||||||||

###Load in data from SpatialPCA drive: https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT
###This does not have ground truth labels
load(HER2_SPCA_preprocessed_filepath) #loads in two objects [location and sp_count]
locations <- data.frame(location)
print(dim(rawcount)) # The count matrix (10,053genes x 607spots) 
print(dim(location)) # The location matrix (607spots x 2)

###Load in raw data and ground truth labels from HER2 github: https://github.com/almaan/her2st/tree/master/data/ST-spotfiles
#Raw data (will be used to index ground truth labels)
H1_raw <- read.table(HER2_H1_raw_filepath, header=TRUE, sep='\t')
#Ground truth labels
truth_df <- read.table(HER2_gt_labels, header=TRUE, sep="\t")





## ||||||||||||||||||||
## Extract ground truth labels ----
## ||||||||||||||||||||

###Reformat ground truth labels
truth_df$x <- round(truth_df$x)
truth_df$y <- round(truth_df$y)

### Filter H1 raw based on colnames of rawcount
H1_raw_filtered <- H1_raw[match(colnames(rawcount), paste0(H1_raw$x, "x", H1_raw$y, "_1")), ]

###Merge (match) ground truth with filtered spots to get final spots with truth labels
truth_df_filtered <- merge(truth_df, H1_raw_filtered, by=c("x", "y"))[c('x','y','label')]

###Add original order column to keep track of location order
locations$original_order <- seq(nrow(locations))

###Merge truth_df_filtered with locations to reorder back to original order
truth_df_filtered <- merge(truth_df_filtered, locations, by=c("x", "y"))
truth_df_filtered <- truth_df_filtered[order(truth_df_filtered$original_order), ]
rownames(truth_df_filtered) <- paste0(truth_df_filtered$x, "x", truth_df_filtered$y, "_1")

###Extract labels, mapping and locations
#labels
truth_lbls <- as.integer(factor(truth_df_filtered$label))
truth_lbls <- as.factor(truth_lbls)
#mapping
mapping <- data.frame(int=levels(truth_lbls), str=levels(factor(truth_df$label)))
#locations
locations_filtered <- truth_df_filtered[,c("x", "y")]
locations_filtered_mat <- as.matrix(locations_filtered)





## ||||||||||||||||||||
## Run Spatial PCA ----
## ||||||||||||||||||||

## Define Time storage variables
trial_list <- list()

tic("Overall time")
###Create SpatialPCA object
tic("Create SPCA object")
SPCA_obj = CreateSpatialPCAObject(counts=rawcount, location=locations_filtered_mat, project = "SpatialPCA",
                            gene.type="spatial",sparkversion="spark", gene.number=3000,
                            customGenelist=NULL,min.loctions = 20, min.features=20)
toc() #Estimate time: 30.5 sec (mac M1 chip)

###Run rest of SpatialPCA pipeline
mem <- peakRAM({
  tic("Building Kernel") 
  SPCA_obj = SpatialPCA_buildKernel(SPCA_obj, kerneltype="gaussian", bandwidthtype="SJ")
  toc() #Estimate time: 0.1 sec (mac M1 chip)
  tic("Estimate Loading")
  SPCA_obj = SpatialPCA_EstimateLoading(SPCA_obj,fast=FALSE,SpatialPCnum=20)
  toc() #Estimate time: 9.5 sec (mac M1 chip)
  tic("Calculating Spatial PCs")
  SPCA_obj = SpatialPCA_SpatialPCs(SPCA_obj, fast=FALSE)
  toc() #Estimate time: 0.5 sec (mac M1 chip)
})
toc() #Estimate time: 41 sec (mac M1 chip)





## ||||||||||||||||||||
## TIME TRIALS: Run Spatial PCA ----
## ||||||||||||||||||||

if (RUN_time_trials){
  n_trials <- 4
  trial_list <- list()
  
  for (i in 1:n_trials){
    ### Clear timing log and print iteration number
    tic.clearlog()
    cat("Iteration:", i, "\n")
    
    tic("Overall time")
    ###Create SpatialPCA object
    tic("Create SPCA object")
    SPCA_obj_timetrial = CreateSpatialPCAObject(counts=rawcount, location=locations_filtered_mat, project = "SpatialPCA",
                                      gene.type="spatial",sparkversion="spark", gene.number=3000,
                                      customGenelist=NULL,min.loctions = 20, min.features=20)
    toc(log=T) #Estimate time: 30.5 sec (mac M1 chip)
    
    ###Run rest of SpatialPCA pipeline
    mem <- peakRAM({
      tic("Building Kernel") 
      SPCA_obj_timetrial = SpatialPCA_buildKernel(SPCA_obj_timetrial, kerneltype="gaussian", bandwidthtype="SJ")
      toc(log=T) #Estimate time: 0.1 sec (mac M1 chip)
      tic("Estimate Loading")
      SPCA_obj_timetrial = SpatialPCA_EstimateLoading(SPCA_obj_timetrial,fast=FALSE,SpatialPCnum=20)
      toc(log=T) #Estimate time: 9.5 sec (mac M1 chip)
      tic("Calculating Spatial PCs")
      SPCA_obj_timetrial = SpatialPCA_SpatialPCs(SPCA_obj_timetrial, fast=FALSE)
      toc(log=T) #Estimate time: 0.5 sec (mac M1 chip)
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
    
  }
  
  SPCA_timetrials_df <- bind_rows(trial_list)
  
  ##Save time trials
  filepath <- glue("{results_dir}SPCA_timetrialresults_{name.dataset}.RData")
  save(SPCA_timetrials_df,
       file=filepath)
  
}








## ||||||||||||||||||||
## Flip coordinates ----
## ||||||||||||||||||||

###Readjust locations because they are flipped across the y axis
locations_final <- locations_filtered
locations_final$y <-  max(locations_filtered$y) + min(locations_filtered$y) - locations_filtered$y 





## ||||||||||||||||||||
## Walktrap CLustering Algorithm ----
## ||||||||||||||||||||

###Run clustering 
clusternum <- 6
knn <- round(sqrt(dim(SPCA_obj@location)[1]))
SpatialPCs <- as.matrix(SPCA_obj@SpatialPCs)
SPCA_clusterlabels <- walktrap_clustering(clusternum = clusternum, latent_dat=SpatialPCs, knearest = knn)
clusterlabel_refine=refine_cluster_10x(SPCA_clusterlabels,SPCA_obj@location,shape="square")





## ||||||||||||||||||||
## Redo labels based off ARI mapping and save ----
## ||||||||||||||||||||

###Align labels
aligned_SPCA_labels <- get_aligned_labels(clustering_labels=clusterlabel_refine, ground_truth=truth_lbls)





## ||||||||||||||||||||
## Calculate Adjusted Rand Index and CHAOS ----
## ||||||||||||||||||||

###Remove "undetermined labels" (label 7) to compute accuracy (only neccessary if clusternum = 7)
idx_remove <- which(truth_lbls == '7')
truth_lbls_no7 <- truth_lbls[-idx_remove]
aligned_SPCA_labels_no7 <- aligned_SPCA_labels[-idx_remove]

###Calculate ARI and CHAOS
SPCA_ARI <- adjustedRandIndex(aligned_SPCA_labels_no7, truth_lbls_no7)
SPCA_CHAOS <- fx_CHAOS(aligned_SPCA_labels_no7, locations_final)





## ||||||||||||||||||||
## Plot spatial components ----
## ||||||||||||||||||||

###Define colors
cbp_spatialpca = c("mediumaquamarine", "chocolate1","dodgerblue",  "#F0E442","palegreen4","lightblue2","plum1")

###Ground truth
plot_cluster(legend="right", location=as.matrix(locations_final), truth_lbls, pointsize=5,
                     text_size=20 ,title_in=paste0("Ground Truth"),color_in=cbp_spatialpca)

###SpatialPCA
plot_cluster(legend="right", location=as.matrix(locations_final), aligned_SPCA_labels, pointsize=5,
                     text_size=20 ,title_in=glue("SpatialPCA (ARI:{round(SPCA_ARI,3)}, CHAOS: {round(SPCA_CHAOS,3)})"), color_in=cbp_spatialpca)





## ||||||||||||||||||||
## Save (or load) files ----
## ||||||||||||||||||||

###Save all objects in one Rdata file
name.dataset <- "HER2"
locations_final_SPCA <- locations_final
filepath <- glue("{results_dir}SPCA_results_{name.dataset}.RData")
save(SPCA_obj, locations_final_SPCA, SPCA_ARI, SPCA_CHAOS, 
     aligned_SPCA_labels, truth_lbls, truth_lbls_no7, 
     file=filepath)










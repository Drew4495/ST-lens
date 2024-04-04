# -------------------------------
# Define filepaths for loading in and saving data
# -------------------------------
###Directories for loading/saving in files/data (REPLACE WITH PERSONAL PATHS)
source_dir <- "~/Codebook/Projects/spatial_transcriptomics_GRSVD/Paper/Manuscript_scripts/src/"
data_dir <- "/Users/aburns/Codebook/Projects/spatial_transcriptomics_GRSVD/Paper/Manuscript_data/Cerebellum/"
cerebellum_data_filepath <- "/Users/aburns/Codebook/Projects/spatial_transcriptomics_GRSVD_R/data/Preprocessed/Mouse_cerebellum_slideseq/mouse_cerebellum_slideseq_data_lulushang.rds"
pseudo_ground_truth_labels_filepath <- "/Users/aburns/Codebook/Projects/spatial_transcriptomics_GRSVD_R/results/Pietro/SlideseqCerebellum_results2/Cerebellum_PseudoTruth_lbls_lambd1eneg3_knn220_npc20.csv"
pseudotruth_coords_filepath <- "/Users/aburns/Codebook/Projects/spatial_transcriptomics_GRSVD_R/results/Pietro/SlideseqCerebellum_results2/Cerebellum_PostFPCA_locations_PIETRO.csv"

# -------------------------------
# Install and load in libraries
# -------------------------------
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

###Load in functions/variables form source (src) files
source(glue("{source_dir}src_support_fxns.R"))



# -------------------------------
# Load in data
# -------------------------------
load(cerebellum_data_filepath) #loads in two objects [location and sp_count]
print(dim(sp_count)) # The count matrix (17,729genes x 25,551spots) 
print(dim(location)) # The location matrix (25,551spots x 2)

###Load pseudotruth labels and coords for alignment
pseudotruth_lbls <- read.csv(pseudo_ground_truth_labels_filepath, header=TRUE)
pseudotruth_coords <- read.csv(pseudotruth_coords_filepath, header=TRUE)



# -------------------------------
# Run Spatial PCA
# -------------------------------
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
  SPCA_obj = SpatialPCA_buildKernel(slideseq, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
  toc()
  tic("Estimate Loading")
  SPCA_obj = SpatialPCA_EstimateLoading(slideseq,fast=TRUE,SpatialPCnum=20)
  toc()
  tic("Calculating Spatial PCs")
  SPCA_obj = SpatialPCA_SpatialPCs(slideseq, fast=TRUE)
  toc()
})
toc()

###Save new locations
locations_postSPCA <- SPCA_obj@location



# -------------------------------
# Louvain CLustering Algorithm
# -------------------------------
###Run clustering
clusternum <- 8
knn <- round(sqrt(dim(SPCA_obj@SpatialPCs)[2]))
SpatialPCs <- as.matrix(SPCA_obj@SpatialPCs)
SPCA_clusterlabels <- louvain_clustering(clusternum = clusternum, latent_dat=SpatialPCs, knearest = knn)

###With 9 clustering
clusternum <- 9
knn <- round(sqrt(dim(SPCA_obj@SpatialPCs)[2]))
SpatialPCs <- as.matrix(SPCA_obj@SpatialPCs)
SPCA_9clusterlabels <- louvain_clustering(clusternum = clusternum, latent_dat=SpatialPCs, knearest = knn)



# -------------------------------
# Align labels based on pseudo ground truth
# -------------------------------
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
aligned_SPCA_labels <- get_aligned_labels(clustering_labels=SPCA_clusterlabels, ground_truth=pseudotruth_lbls_matched)



# -------------------------------
# Calculate PAS and CHAOS scores and save
# -------------------------------
###Calculate CHAOS and PAS
CHAOS_score <- fx_CHAOS(aligned_SPCA_labels, locations_postSPCA)
PAS_score <- fx_PAS(aligned_SPCA_labels, locations_postSPCA)
CHAOS_PAS_df <- data.frame(CHAOS=CHAOS_score, PAS=PAS_score)



# -------------------------------
# Plot spatial components
# -------------------------------
###Define colors
cbp=c("#F5F5DC","#FF7070" , "#C3FFFF" ,"#00E1A1","#256BAA","#FED439" ,"#4DBBD5", "#87CEEB")

###Pseudo-Ground_truth
plot_cluster(location=pseudotruth_coords, as.character(pseudotruth_lbls[,1]), pointsize=1.5,
             title_in="Pseudotruth Labels",color_in=cbp,legend="right")

###SpatialPCA
plot_cluster(location=locations_postSPCA, aligned_SPCA_labels, pointsize=1.5,
             title_in=glue("SpatialPCA (CHAOS: {round(CHAOS_score,3)}, PAS: {round(PAS_score,3)})"),color_in=cbp,legend="right")



# -------------------------------
# Save/Load files
# -------------------------------
###Save SpatialPCA object
filepath <- glue("{data_dir}Cerebellum_SPCA_obj.Rdata")
#save(SPCA_obj, file=filepath)
load(file = filepath)

###Save normalized expression matrix
filepath = glue("{data_dir}Cerebellum_gene_by_spot_mat_normalized_zero_mean_1stddev.csv")
#write.csv(SPCA_obj@normalized_expr, file=filepath, row.names=TRUE)

###Save raw xy coords
filepath = glue("{data_dir}SpatialPCA_Cerebellum_xy_coords_raw.csv")
#write.csv(location, file=filepath, row.names=TRUE)
locations_final <- read.csv(file = filepath, header = TRUE, row.names=1)

###Save Post SPCA xy coords
filepath = glue("{data_dir}SpatialPCA_Cerebellum_xy_coords_postSPCA.csv")
#write.csv(locations_postSPCA, file=filepath, row.names=TRUE)
locations_final <- read.csv(file = filepath, header = TRUE, row.names=1)

###Save spatial components from HER2 object
filepath = "{data_dir}SPCA_Cerebellum_SpatialComps.csv"
#write.csv(SPCA_obj@SpatialPCs, file=filepath, row.names=TRUE)

###Save ALIGNED cluster labels for SPCA
filepath = glue("{data_dir}SPCA_Cerebellum_ALIGNEDclusterlabels_wolftrap.csv")
#write.csv(aligned_SPCA_labels, file=filepath, row.names=FALSE)
aligned_SPCA_labels <- as.character(read.csv(file=filepath, header=TRUE)[,1])

####Save ARI and CHAOS score
CHAOS_PAS_df <- data.frame(CHAOS=CHAOS_score, PAS=PAS_score)
filepath <- glue("{data_dir}SpatialPCA_Cerebellum_CHAOSandPAS_score.csv")
#write.csv(CHAOS_PAS_df, filepath, row.names = FALSE)
CHAOS_PAS_df <- read.csv(file=filepath)

###Save all objects in one Rdata file
name.dataset <- "SlideseqCerebellum"
locations_final_SPCA <- locations_final
SPCA_CHAOS <- CHAOS_score
SPCA_PAS <- PAS_score
filepath <- glue("{data_dir}SPCA_results_{name.dataset}.RData")
save(SPCA_obj, locations_final_SPCA, SPCA_PAS, SPCA_CHAOS, aligned_SPCA_labels, 
     SPCA_9clusterlabels, file=filepath)





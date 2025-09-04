## Workflow for manuscript

Order: **1a) Meshing → 1b) Preprocessing → 2) Analysis (Functional PCA) → 3) Decide # Components → 4) Clustering → 5) Visualization**

Meshing and Preprocessing scripts are interchangeable. Both need to be completed before "Analysis" script is run.
For specific scripts used in HER2, DLPFC, and Cerebellum datasets, in the corresponding manuscript, consult the """" in """"


## Script Descriptions (Inputs & Outputs)

Each dataset has the corresponding scripts

- **Meshing**  
  Build a hexagonal lattice over the tissue, simplify the boundary, and generate a **finite-element (FEM)** triangular mesh that matches the tissue shape.  
  _Output:_ `data/{name_dataset}/mesh.RData` (contains `mesh` and `lattice`)  
  _Figures:_ step-by-step PDFs in `images/DLPFC/meshing/` named `0_...pdf` through `9_...pdf`

- **Preprocessing**  
  Use SpatialPCA’s preprocessing to select spatially informative genes, filter spots/genes, normalize counts, then **z-score each gene** (row-wise `scale()`).  
  _Output:_ `data/{name_dataset}/preprocessed_data.RData`  
  _(Includes both the original “snapshots” and the processed `locations`, z-scored `counts`, names, and `true_labels`.)_

- **Analysis (Functional PCA = fPCA)**  
  Align data to the mesh domain (drop out-of-domain spots) and fit **fPCA** with **generalized cross-validation (GCV)** to choose smoothing. Save mean and component fields and evaluate them on a high-resolution grid and at spot locations.  
  _Outputs:_  
  - `results/{name_dataset}/fPCA.RData` (mean, scores, loadings, evaluated fields, timing)  
  - `data/{name_dataset}/analyzed_data.RData` (domain-aligned data; not fPCA results)  
  - _Figures:_ `images/DLPFC/analysis/1_HR_grid.pdf`, `2_mean.pdf`, `2_mean_at_locations.pdf`, `3_components_at_locations.pdf`  
  - _Note for timing scripts:_ `results/{name_dataset}/domain_input_for_FPCA.RData` (contains `mesh`, `counts`, `locations`; move the `save()` after grid creation if you want `grid` included)

- **Decide number of spatial components**  
  Plot the **normalized reconstruction error** (RMSE of data vs. fPCA low-rank reconstructions) across component counts and pick the elbow.  
  _Outputs:_  
  - `images/{name_dataset}/NSC_optimization/Optimal_nComp_Selection.pdf`  
  - `results/{name_dataset}/nComp_opt.RData` (e.g., `nComp_opt <- 3`)  

- **Clustering**  
  Run **Walktrap community detection** on the first `nComp_opt` fPCA components on a high-resolution grid, then assign each original spot the label of its nearest grid point. Score with **Adjusted Rand Index (ARI)** vs ground truth and **align cluster IDs** to ground-truth IDs via a greedy ARI-based mapping. (Some may repeat for differnts NSC numbers for comparison.)  
  _Outputs (per NSC):_  
  - `results/{name_dataset}/fPCA_clustering_results_GCV_NoMean_{NSC}.RData` (lists of HR/LR labels and ARIs, keyed by `knearest_{K}`)  
  - `results/{name_dataset}/fPCA_clustering_results_ALIGNED_GCV_NoMean_{NSC}.RData` (aligned versions)  
  - _Figures:_ per-KNN PDFs in `images/DLPFC/clustering/GCV_NoMean/...`:  
    `HR_clustering_k{K}.pdf`, `clustering_k{K}.pdf`, plus `aggregated_KNN_HR_plot.pdf` and `aggregated_KNN_LR_plot.pdf`  

- **Visualization**  
  Export manuscript-style figures: **Ground Truth vs SpatialPCA vs fPCA** label maps (with ARI/CHAOS annotations) and selected spatial component maps (scaled and unscaled).  
  _Outputs:_ PDFs in `images/{name_dataset}/manuscript/` (e.g., `DLPFC_sample9_GroundTruth.pdf`, `DLPFC_sample9_SPCA.pdf`, `DLPFC_sample9_FPCA.pdf`, and component grids)

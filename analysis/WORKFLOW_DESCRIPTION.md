# Workflow for manuscript

Order: **1a) Meshing → 1b) Preprocessing → 2) Analysis (Functional PCA) → 3) Decide # Components → 4) Clustering → 5) Visualization**

Meshing and Preprocessing scripts are interchangeable. Both need to be completed before "Analysis" script is run.


## Script Descriptions (Inputs & Outputs)

Each dataset has the corresponding scripts

- ### `Meshing_{name_dataset}.R`  
  Build a hexagonal lattice over the tissue, simplify the boundary, and generate a **finite-element (FEM)** triangular mesh that matches the tissue shape.

  **Inputs**
  - `data/{name_dataset}/data.RData` (expects `counts`, `locations`, `true_labels`)

  **Output**
  - `data/{name_dataset}/mesh.RData` (contains `mesh` and `lattice`)

  **Figures**
  - Step-by-step PDFs in `images/{name_dataset}/meshing/` named `0_...pdf` through `9_...pdf`

- ### `Preprocessing_{name_dataset}.R` 
  preprocessing to select spatially informative genes, filter spots/genes, normalize counts, then **z-score each gene** (row-wise `scale()`).

  **Inputs**
  - `data/{name_dataset}/data.RData`

  **Output**
  - `data/{name_dataset}/preprocessed_data.RData`  
    (Includes both the original snapshots and the processed `locations`, z-scored `counts`, names, and `true_labels`.)

- ### `Analysis_{name_dataset}.R` (Functional PCA = fPCA) 
  Align data to the mesh domain (drop out-of-domain spots) and fit **fPCA** with **generalized cross-validation (GCV)** to choose smoothing. Save mean and component fields and evaluate them on a high-resolution grid and at spot locations.

  **Inputs**
  - `data/{name_dataset}/preprocessed_data.RData`
  - `data/{name_dataset}/mesh.RData`

  **Outputs**
  - `results/{name_dataset}/fPCA.RData` (mean, scores, loadings, evaluated fields, timing)
  - `data/{name_dataset}/analyzed_data.RData` (domain-aligned data; not fPCA results)

  **Figures**
  - `images/{name_dataset}/analysis/1_HR_grid.pdf`
  - `images/{name_dataset}/analysis/2_mean.pdf`
  - `images/{name_dataset}/analysis/2_mean_at_locations.pdf`
  - `images/{name_dataset}/analysis/3_components_at_locations.pdf`


- ### `decide_NSC_opt_{name_dataset}.R` (Decide number of spatial components)  
  Plot the **normalized reconstruction error** (RMSE of data vs. fPCA low-rank reconstructions) across component counts and pick the elbow.

  **Inputs**
  - `data/{name_dataset}/mesh.RData`
  - `data/{name_dataset}/analyzed_data.RData`
  - `results/{name_dataset}/fPCA.RData`

  **Outputs**
  - `images/{name_dataset}/NSC_optimization/Optimal_nComp_Selection.pdf`
  - `results/{name_dataset}/nComp_opt.RData` (e.g., `nComp_opt <- 3`)

- ### `clustering_{name_dataset}.R` 
  Run **Walktrap community detection** on the first `nComp_opt` fPCA components on a high-resolution grid, then assign each original spot the label of its nearest grid point. Score with **Adjusted Rand Index (ARI)** vs ground truth and **align cluster IDs** to ground-truth IDs via a greedy ARI-based mapping. (Optionally repeat for other NSC values, e.g., 8 and 2, for comparison.)

  **Inputs**
  - `data/{name_dataset}/mesh.RData`
  - `data/{name_dataset}/analyzed_data.RData`
  - `results/{name_dataset}/fPCA.RData`
  - `results/{name_dataset}/nComp_opt.RData`

  **Outputs (per NSC)**
  - `results/{name_dataset}/fPCA_clustering_results_GCV_NoMean_{NSC}.RData` (lists of HR/LR labels and ARIs, keyed by `knearest_{K}`)
  - `results/{name_dataset}/fPCA_clustering_results_ALIGNED_GCV_NoMean_{NSC}.RData` (aligned versions)

  **Figures**
  - Per-KNN PDFs in `images/{name_dataset}/clustering/GCV_NoMean/...`:  
    `HR_clustering_k{K}.pdf`, `clustering_k{K}.pdf`, plus `aggregated_KNN_HR_plot.pdf` and `aggregated_KNN_LR_plot.pdf`

- ### `SPCA_and_FPCA_visualization_{name_dataset}.R`  
  Export manuscript-style figures: **Ground Truth vs SpatialPCA vs fPCA** label maps (with ARI/CHAOS annotations) and selected spatial component maps (scaled and unscaled).

  **Inputs**
  - `data/{name_dataset}/data.RData`
  - `results/{name_dataset}/SPCA_allresults_NEW_{name_dataset}.RData`
  - `data/{name_dataset}/analyzed_data.RData`
  - `results/{name_dataset}/fPCA.RData`
  - `results/{name_dataset}/fPCA_clustering_results_ALIGNED_GCV_NoMean_3.RData`

  **Outputs**
  - PDFs in `images/{name_dataset}/manuscript/` (e.g., `{name_dataset}_GroundTruth.pdf`, `{name_dataset}_SPCA.pdf`, `{name_dataset}_FPCA.pdf`, and component grids)





# Accessory Scripts

- ### `fpca_timetrials.R`
  - **What it does:** Benchmarks runtime for (1) meshing from locations and (2) fPCA on a prepackaged mesh+data bundle. Saves tidy timing tables.
  - **Inputs:**
    - `data/{name_dataset}/data.RData` — for meshing trials (uses `locations`)
    - `results/{name_dataset}/domain_input_for_FPCA.RData` — for fPCA trials (contains `mesh`, `counts`, `locations`)
  - **Outputs:**
    - `results/{name_dataset}/timetrials_mesh_df.RData` — data frame `mesh_timetrials_df` (per-trial elapsed times)
    - `results/{name_dataset}/timetrials_fpca_df.RData` — data frame `fpca_timetrials_df` (per-trial elapsed times)
  - **Notes:**
    - Meshing trials **do not** save a mesh; they only save timing results.

---

- ### `SpatialPCA.R`
  - **What it does:** Runs SpatialPCA preprocessing + model, optionally repeats for timing, then clusters SpatialPCA components, refines labels, computes ARI and CHAOS, and saves a compact results bundle.
  - **Inputs:**
    - `data/{name_dataset}/data.RData` — raw data (`counts`, `locations`, `true_labels`)
    - helper functions from `src/`
  - **Outputs:**
    - `results/{name_dataset}/SPCA_allresults_NEW_{name_dataset}.RData`  
      *(Contains `SPCA_obj`, `xy_coords` (raw), `locations_final_SPCA` (filtered), `aligned_SPCA_labels`, `SPCA_ARI`, `SPCA_CHAOS`.)*
    - *(Optional, if `RUN_time_trials <- TRUE`)* `results/{name_dataset}/SPCA_timetrialresults_{name_dataset}.RData` — timing data frame `SPCA_timetrials_df`



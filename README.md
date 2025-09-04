# fPCA for Spatial Transcriptomics

Functional Principal Component Analysis (fPCA) is a fast, accurate workflow for spatial domain detection in spatial transcriptomics that models gene expression as smooth functions over a finite-element mesh—preserving complex tissue geometry (holes, non-convexities) and true spatial relationships. 
 
Compared with leading methods, fPCA runs substantially faster, achieves higher accuracy which was tested with simulated datasets with known ground truth.

It also provides inherent high-resolution reconstructions and scales to dense, high-coverage assays; results are demonstrated on HER2+ breast cancer, human DLPFC, and mouse cerebellum. 
 

![fPCA Overview](README_images/ST_Method_Overview_Vertical.png)

A non-convex shape is one that in which at least there exists a pair of points in which a straight line segment between them is contained completely within the shape itself. In the example below, the red lines satisfy this property while the blue line does not, indicating the right shape is non-convex.
![Convexity and Non-Convexity Example](README_images/convex.png)

Photo adapted from https://undergroundmathematics.org/glossary/convex-shape




## Tutorials

Tutorials created by Marco Galliani (acknowledged in manuscript)
- [fPCA 2D Tutorial (code on GitHub)](https://github.com/Drew4495/ST-lens/blob/main/Tutorials/fPCA_2D.html)  
- [fPCA 2D Tutorial (rendered)](https://htmlpreview.github.io/?https://raw.githubusercontent.com/Drew4495/ST-lens/main/Tutorials/fPCA_2D.html)




## Pipeline (at a glance)

Order: **Meshing → Preprocessing → Analysis (Functional PCA) → Decide # Components → Clustering → Visualization**

- **Meshing**  
  Build a hexagonal lattice over the tissue, simplify the boundary, and generate a **finite-element (FEM)** triangular mesh that matches the tissue shape.  
  _Output:_ `data/DLPFC/mesh.RData` (contains `mesh` and `lattice`)  
  _Figures:_ step-by-step PDFs in `images/DLPFC/meshing/` named `0_...pdf` through `9_...pdf`

- **Preprocessing**  
  Use SpatialPCA’s preprocessing to select spatially informative genes, filter spots/genes, normalize counts, then **z-score each gene** (row-wise `scale()`).  
  _Output:_ `data/DLPFC/preprocessed_data.RData`  
  _(Includes both the original “snapshots” and the processed `locations`, z-scored `counts`, names, and `true_labels`.)_

- **Analysis (Functional PCA = fPCA)**  
  Align data to the mesh domain (drop out-of-domain spots) and fit **fPCA** with **generalized cross-validation (GCV)** to choose smoothing. Save mean and component fields and evaluate them on a high-resolution grid and at spot locations.  
  _Outputs:_  
  - `results/DLPFC/fPCA.RData` (mean, scores, loadings, evaluated fields, timing)  
  - `data/DLPFC/analyzed_data.RData` (domain-aligned data; not fPCA results)  
  - _Figures:_ `images/DLPFC/analysis/1_HR_grid.pdf`, `2_mean.pdf`, `2_mean_at_locations.pdf`, `3_components_at_locations.pdf`  
  - _Note for timing scripts:_ `results/DLPFC/domain_input_for_FPCA.RData` (contains `mesh`, `counts`, `locations`; move the `save()` after grid creation if you want `grid` included)

- **Decide number of spatial components**  
  Plot the **normalized reconstruction error** (RMSE of data vs. fPCA low-rank reconstructions) across component counts and pick the elbow.  
  _Outputs:_  
  - `images/DLPFC/NSC_optimization/Optimal_nComp_Selection.pdf`  
  - `results/DLPFC/nComp_opt.RData` (e.g., `nComp_opt <- 3`)  

- **Clustering**  
  Run **Walktrap community detection** on the first `nComp_opt` fPCA components on a high-resolution grid, then assign each original spot the label of its nearest grid point. Score with **Adjusted Rand Index (ARI)** vs ground truth and **align cluster IDs** to ground-truth IDs via a greedy ARI-based mapping. (Also repeats for NSC = 8 and 2 for comparison.)  
  _Outputs (per NSC):_  
  - `results/DLPFC/fPCA_clustering_results_GCV_NoMean_{NSC}.RData` (lists of HR/LR labels and ARIs, keyed by `knearest_{K}`)  
  - `results/DLPFC/fPCA_clustering_results_ALIGNED_GCV_NoMean_{NSC}.RData` (aligned versions)  
  - _Figures:_ per-KNN PDFs in `images/DLPFC/clustering/GCV_NoMean/...`:  
    `HR_clustering_k{K}.pdf`, `clustering_k{K}.pdf`, plus `aggregated_KNN_HR_plot.pdf` and `aggregated_KNN_LR_plot.pdf`  

- **Visualization**  
  Export manuscript-style figures: **Ground Truth vs SpatialPCA vs fPCA** label maps (with ARI/CHAOS annotations) and selected spatial component maps (scaled and unscaled).  
  _Outputs:_ PDFs in `images/DLPFC/manuscript/` (e.g., `DLPFC_sample9_GroundTruth.pdf`, `DLPFC_sample9_SPCA.pdf`, `DLPFC_sample9_FPCA.pdf`, and component grids)

👉 For full inputs/outputs per script, folder structure, and caveats, see **[Full Workflow & File Map](docs/WORKFLOW.md)**.













TEMPLATE
# Foobar

Foobar is a Python library for dealing with word pluralization.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)

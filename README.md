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




## Generic Pipeline (at a glance)

Order: **1a) Meshing → 1b) Preprocessing → 2) Analysis (Functional PCA) → 3) Decide # Components → 4) Clustering → 5) Visualization**

Meshing and Preprocessing scripts are interchangeable. Both need to be completed before "Analysis" script is run.
For specific scripts used in HER2, DLPFC, and Cerebellum datasets, in the corresponding manuscript, consult [WORKFLOW_DESCRIPTION.md](analysis/WORKFLOW_DESCRIPTION.md) in `analysis/`

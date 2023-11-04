# Optimal components ----
# |||||||||||||||||||||||

# Model
results_fPCA <- fPCA(mesh, locations, counts.normalized, lambdas, nComp)

# Results
scores <- results_fPCA$scores
loadings <- results_fPCA$loadings
loadings_nodes <- results_fPCA$loadings_nodes


## Save components ----
## ||||||||||||||

save(scores,
     loadings, loadings_nodes,
     # Saving options
     file = paste(directory.results, name.dataset, "_components", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(scores, 
   loadings, loadings_nodes,
   results_fPCA)

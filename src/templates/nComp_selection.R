# nComp selection ----
# ||||||||||||||||||||

# Model
results_fPCA <- fPCA_CS(mesh, locations, counts.normalized, lambda, nComp)

# Results
scores <- results_fPCA$scores
loadings <- results_fPCA$loadings
loadings_nodes <- results_fPCA$loadings_nodes

# Residuals
norm <- RMSE(counts.normalized,0)
residuals <- fPCA_residuals(counts.normalized, scores, loadings)
residuals_norm <- list()
residuals_norm[["N"]] <- 1
residuals_norm <- c(residuals_norm, residuals/norm)

## Save results ----
## |||||||||||||||||

save(residuals_norm, nComp_opt, scores, loadings, loadings_nodes,
     # Saving options
     file = paste(directory.results, name.dataset, "_nComp_selection", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(results_fPCA, 
   scores, loadings, loadings_nodes,
   norm, residuals, residuals_norm)

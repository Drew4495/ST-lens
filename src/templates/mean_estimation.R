# Mean estimation ----
# ||||||||||||||||||||

# Model
results_FRPDE <- FRPDE(mesh, locations, counts, lambdas)

# Results
counts.mean <- results_FRPDE$fitted
counts.mean_nodes <- results_FRPDE$fitted_nodes
lambda_opt <- results_FRPDE$lambda_opt

# Data centering
counts.centered <- counts - rep(1, nrow(counts))%*%t(counts.mean)


## Save mean ----
## ||||||||||||||

save(counts.mean, counts.mean_nodes,
     counts.centered,
     lambda_opt,
     # Saving options
     file = paste(directory.results, name.dataset, "_mean", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(counts.mean, counts.mean_nodes,
   counts.centered,
   results_FRPDE, lambda_opt)


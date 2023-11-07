# fdaPDE ----
# |||||||||||

## FEM utilities ----
## ||||||||||||||||||

mesh_data <- function(mesh) {
  
  nodes <- data.frame(mesh$nodes)
  colnames(nodes) <- c("x", "y")
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  return(mesh_data)
}

pde.init <- function(mesh) {
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data(mesh))
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  return(pde)
}

field.eval <- function(grid, f_nodes, mesh) {
  
  FEMbasis = create.FEM.basis(mesh)
  
  FEMfunction = FEM(f_nodes, FEMbasis)
  f_grid <- eval.FEM(FEMfunction, grid)
  
  return(f_grid)
}


## Models ----
## |||||||||||

### Utilities ----
### ||||||||||||||

RMSE <- function(x1, x2){
  sqrt(sum((x1 - x2)^2)/length(x1))
}

fPCA_residuals <- function(X, scores, loadings) {
  
  residuals <- c()
  for(h in 1:ncol(scores)){
    residuals[h] <- RMSE(X, scores[,1:h] %*% t(loadings[,1:h]))
  }
  names(residuals) <- as.character(1:ncol(scores))
  
  return(residuals)
}

fPCA_normalize <- function(scores, loadings, loadings_nodes) {
  
  for(h in 1:ncol(scores)){
    norm <- sqrt(sum(scores[,h]^2))
    loadings[,h] <- loadings[,h]*norm
    loadings_nodes[,h] <- loadings_nodes[,h]*norm
    scores[,h] <- scores[,h]/norm
  }
  
  return(list(scores = scores,
              loadings = loadings,
              loadings_nodes = loadings_nodes))
}


### Wrappers ----
### |||||||||||||

# Functional Regression with PDE regularization model

FRPDE <- function(mesh, locations, X, lambdas, verbose = FALSE){
  
  # Define and init model
  model_FRPDE <- new(FRPDE_Laplacian_2D_GeoStatLocations, pde.init(mesh))
  model_FRPDE$set_verbose(verbose)
  
  # Set locations
  model_FRPDE$set_locations(as.matrix(locations))
  
  # Set observations
  model_FRPDE$set_observations(as.matrix(X))
  
  # Tune
  lambda_opt <- model_FRPDE$tune(lambdas)
  
  # Solve
  model_FRPDE$set_lambda_s(lambda_opt)
  model_FRPDE$solve()
  
  # Results
  fitted <- t(model_FRPDE$fitted())
  fitted_nodes <- t(model_FRPDE$f())
  
  return(list(lambda_opt = lambda_opt,
              fitted = fitted,
              fitted_nodes = fitted_nodes))
}

# Functional Principal Component Analysis (Closed-Form solution implementation)

fPCA_CS <- function(mesh, locations, X, lambda, nComp, verbose = FALSE){
  
  # Define and init model
  model_fPCA <- new(FPCA_CS_Laplacian_2D_GeoStatLocations, pde.init(mesh))
  model_fPCA$set_npc(nComp)
  model_fPCA$set_lambda_s(lambda)
  
  # Model options
  model_fPCA$set_mass_lumping(TRUE)
  model_fPCA$set_coefficients_position(2)
  model_fPCA$set_verbose(verbose)
  
  # Set locations
  model_fPCA$set_locations(as.matrix(locations))
  
  # Set observations
  model_fPCA$set_observations(as.matrix(X))
  
  # Solve
  model_fPCA$solve()
  
  # Results
  scores <- model_fPCA$scores()
  loadings <- model_fPCA$loadings()
  loadings_nodes <- model_fPCA$W()
  
  return(list(scores = scores,
              loadings = loadings,
              loadings_nodes = loadings_nodes))
}

# Functional Principal Component Analysis (Iterative scheme)

fPCA <- function(mesh, locations, X, lambdas, nComp, verbose = FALSE){
  
  # Define and init model
  model_fPCA <- new(FPCA_Laplacian_2D_GeoStatLocations, pde.init(mesh))
  model_fPCA$set_npc(nComp)
  model_fPCA$set_lambdas(lambdas)
  
  # Set locations
  model_fPCA$set_locations(as.matrix(locations))
  
  # Set observations
  model_fPCA$set_observations(as.matrix(X))
  
  # Solve
  model_fPCA$solve()
  
  # Results
  scores <- model_fPCA$scores()
  loadings <-  model_fPCA$loadings()
  loadings_nodes <- model_fPCA$W()
  
  # Normalize results
  results_norm <- fPCA_normalize(scores, loadings, loadings_nodes)
  scores <- results_norm$scores
  loadings <-  results_norm$loadings
  loadings_nodes <- results_norm$loadings_nodes
  
  return(list( scores = scores,
               loadings = loadings,
               loadings_nodes = loadings_nodes))
}
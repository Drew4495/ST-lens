square_grid <- function(bbox, h, seed_point = NULL) {
  
  xmin <- bbox[1,1]
  xmax <- bbox[1,2]
  ymin <- bbox[2,1]
  ymax <- bbox[2,2]
  
  x <- seq(xmin, xmax, by = h)
  y <- seq(ymin, ymax, by = h)
  grid <- expand.grid(x = x, y = y)
  
  return(grid)
}

evaluate_field <- function(grid, f_at_nodes, mesh) {
  
  FEMbasis = create.FEM.basis(mesh)
  
  FEMfunction = FEM(f_at_nodes, FEMbasis)
  f_at_grid <- eval.FEM(FEMfunction, grid)
  
  return(f_at_grid)
}

dist_point_from_points <- function(p, points){
  points <- data.frame(points)
  points$x <- points$x - p$x
  points$y <- points$y - p$y
  dist <- sqrt(points$x^2 + points$y^2)
  
  return(dist)
}

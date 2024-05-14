# plots ----


## mesh generation ----


## plots the discretized domain

plot.discretized_domain <- function(lattice, size = 1) {
  
  # Data
  domain <- lattice$domain        # Discretized domain \Omega_h
  points.b <- lattice$points.b    # Outer boundary points
  points.h <- lattice$points.h    # Holes boundary points
  
  # Colors
  cols.color = c("Outer Boundary" = adjustcolor("red"),
                 "Holes Boundary" = adjustcolor("orange"))
  
  plot <- ggplot() +
    xlab("x") + ylab("y") +
    geom_sf(data = st_as_sf(domain), fill = "grey", alpha = 0.5, color = "black") +
    geom_sf(data = st_as_sf(points.b), aes(color = "Outer Boundary"), size = size) +
    scale_color_manual(name = "", values = cols.color)
  
  # Add holes
  if(pol_holes_number(domain) >= 1){
    for(i in 1:length(points.h)){
      plot <- plot +
        geom_sf(data = st_as_sf(points.h[[i]]), aes(color = "Holes Boundary"), size = size)
    }
  }
  
  return(plot)
}


# plots the mesh:
# - Boundaries are depicted in red
# - Elements are depicted in black

plot.fdaPDE_mesh <- function(mesh) {
  
  ## nodes
  nodes <- data.frame(mesh$nodes)
  colnames(nodes) <- c("x", "y")
  
  ## triangles
  triangles <- data.frame(mesh$triangles)
  colnames(triangles) <- c("x1", "x2", "x3")
  triangles$id <- 1:dim(triangles)[1]
  triangles <- reshape(triangles,
                       idvar = "id",
                       varying = c("x1", "x2", "x3"),
                       v.names = "index",
                       timevar = "order",
                       times = c(1:3),
                       direction = "long")
  triangles_cordinates <- nodes[triangles$index,]
  triangles = cbind(triangles, triangles_cordinates)
  polygons <- list()
  for(id in 1:max(triangles$id)){
    points <- triangles[triangles$id == id, c("x", "y")]
    points <- rbind(points, points[1,])
    polygons[[id]] <- Polygon(points)
  }
  triangles <- SpatialPolygons(list(Polygons(polygons, ID = "elements")))
  
  ## boundary segments
  boundaries <- data.frame(mesh$segments[mesh$segmentsmarkers == 1,])
  boundaries$id <- 1:dim(boundaries)[1]
  boundaries <- reshape(boundaries,
                        idvar = "id",
                        varying = c("X1", "X2"),
                        v.names = "index",
                        timevar = "order",
                        times = c(1:2),
                        direction = "long")
  boundaries_cordinates <- nodes[boundaries$index,]
  boundaries <- cbind(boundaries, boundaries_cordinates)
  lines <- list()
  for(id in 1:max(boundaries$id)){
    points <- boundaries[boundaries$id == id, c("x", "y")]
    lines[[id]] <- Lines(list(Line(points)), ID = id)
  }
  boundaries <- SpatialLines(lines)
  
  ## mesh
  plot <- ggplot() + standard_plot_settings_fields() +
    geom_sf(data = st_as_sf(triangles), color = "black", fill = "transparent",
            linewidth = 0.1) +
    geom_sf(data = st_as_sf(boundaries), color = "red", fill = "transparent")
  
  return(plot)
  
}

## plots the final location on the domain highlighting the discarded ones in red

plot.final_locations <- function(locations, locations.final, lattice, size = 1) {
  
  ## colors
  cols.color <- c("Final" = "black", "Discarded" = "red")
  
  plot <- ggplot() +
    standard_plot_settings_fields() + 
    geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black") +
    geom_sf(data = st_as_sf(locations), aes(color = "Discarded"), size = size) +
    geom_sf(data = st_as_sf(locations.final), aes(color = "Final"), size = size) +
    scale_color_manual(name = "", values = cols.color)
  
  ## remove the legend if there are no discarded points
  if(nrow(locations.final@coords) < nrow(locations@coords)){
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}
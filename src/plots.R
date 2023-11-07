# Plots ----
# ||||||||||

## Settings ----
## |||||||||||||

mesh_generation_examples_plot_settings <- function(){
  theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank())
  
}


## Mesh generation ----
## ||||||||||||||||||||

# Plots the discretized domain

plot.discretized_domain <- function(lattice, size = 1) {
  
  # Data
  domain <- lattice$domain        # Discretized domain \Omega_h
  points.b <- lattice$points.b    # Outer boundary points
  points.h <- lattice$points.h    # Holes boundary points
  
  # Colors
  cols.color = c("Outer Boundary" = adjustcolor("red"),
                 "Holes Boundary" = adjustcolor("orange"))
  
  plot <- ggplot() + standard_plot_settings() +
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

# Plots the mesh:
# - Boundaries are depicted in red
# - Elements are depicted in black

plot.fdaPDE_mesh <- function(mesh) {
  
  # Nodes
  nodes <- data.frame(mesh$nodes)
  colnames(nodes) <- c("x", "y")
  
  # Triangles
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
  
  # Boundary segments
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
  
  # Mesh
  plot <- ggplot() + standard_plot_settings() +
    xlab("x") + ylab("y") +
    geom_sf(data = st_as_sf(triangles), color = "black", fill = "transparent",
            linewidth = 0.1) +
    geom_sf(data = st_as_sf(boundaries), color = "red", fill = "transparent")
  
  return(plot)
  
}

# Plots the final location on the domain highlighting the discarded ones in red

plot.final_locations <- function(locations, locations.final, lattice, size = 1) {
  
  # Colors
  cols.color <- c("Final" = "black", "Discarded" = "red")
  
  plot <- ggplot() +
    xlab("x") + ylab("y") +
    standard_plot_settings() + 
    geom_sf(data = st_as_sf(lattice$domain), fill = "grey", alpha = 0.5, color = "black") +
    geom_sf(data = st_as_sf(locations), aes(color = "Discarded"), size = size) +
    geom_sf(data = st_as_sf(locations.final), aes(color = "Final"), size = size) +
    scale_color_manual(name = "", values = cols.color)
  
  # Remove the legend if there are no discarded points
  if(nrow(locations.final@coords) < nrow(locations@coords)){
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}


## Analysis -----
##|||||||||||||||

# Plots the functional components

plot.components <- function(locations, loadings, type = "points",
                            size = 1,  colormap = "D", limits = NULL,
                            ncol = 5) {
  
  plots <- list()
  for(i in 1:ncol(loadings)){
    if(type == "points"){
      plots[[i]] <- plot.field_points(locations, loadings[,i], size = size,
                                      colormap = colormap, limits = limits) 
      plots[[i]] <- plots[[i]] + guides(color = "none")
    } else if(type == "tile"){
      plots[[i]] <- plot.field_tile(locations, loadings[,i],
                                    colormap = colormap, limits = limits)
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
    plots[[i]] <- plots[[i]] + ggtitle(paste("Comp", i)) + xlab("") + ylab("")
  }
  
  plot <- arrangeGrob(grobs = plots, ncol = min(ncol, nComp))
  
  return(plot)
}

# Guide plot to decide the optimal number of components

plot.nComp_selection <- function(residuals_norm, nComp, nComp_opt) {
  
  # Percentage of data explained
  data_plot <- data.frame(id = 0:nComp,
                          name = factor(names(residuals_norm), levels = names(residuals_norm)),
                          residuals_norm = unlist(residuals_norm))
  plot1 <- ggplot(data = data_plot) +
    geom_hline(yintercept = min(data_plot$residuals_norm), linetype = "dashed", color = "red", size = 0.5) +
    geom_line(aes(x = id, y = residuals_norm, group = 1)) +
    geom_point(aes(x = id, y = residuals_norm), size = 2) +
    geom_point(aes(x = id[1+nComp_opt], y = residuals_norm[1+nComp_opt]), size = 2, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Percentage explained") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings()
  
  # Information gained
  data_plot <- data.frame(id = 0:nComp,
                          name = factor(names(residuals_norm), levels = names(residuals_norm)),
                          gain = - unlist(residuals_norm) + c(NaN, unlist(residuals_norm)[1:(nComp)]))
  plot2 <- ggplot(data = data_plot) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    geom_line(aes(x = id, y = gain, group = 1)) +
    geom_point(aes(x = id, y = gain), size = 2) +
    geom_point(aes(x = id[1+nComp_opt], y = gain[1+nComp_opt]), size = 2, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Percentage gain") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings()
  
  # Title
  title <- textGrob("Optimal nComp selection",
                    gp = gpar(fontsize = 18, fontface = 'bold'))
  
  # Plot
  plot <- arrangeGrob(title, plot1, plot2, nrow = 3, heights = c(1, 4, 4))
  
  return(plot)
}

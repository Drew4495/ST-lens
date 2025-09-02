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








## |||||||||||||||||||||||||||||||||||||||      
# Plotting Functions by Drew Burns based off Pietro Donelli's Code ----
## ||||||||||||||||||||||||||||||||||||||| 


## |||||||||||||
## Settings ----
## |||||||||||||

standard_plot_settings <- function(){
  #' Standard plot settings for ggplot2
  #'
  #' This function applies a minimal theme with customized text settings and legend position for ggplot2 plots.
  #' It is designed to provide a base setting for creating clean and professional-looking plots.
  #'
  #' @return A ggplot2 theme object with specified plot settings.
  #' @import ggplot2
  standard_plot_settings <- theme_minimal() +
    theme(text = element_text(size = 12),
          plot.title = element_text(
            color = "black",
            face = "bold",
            size = 14,
            hjust = 0.5,
            vjust = 1),
          legend.position = "bottom",
          legend.title = element_blank())
}

manuscript_plot_settings <- function(){
  #' Manuscript plot settings for ggplot2
  #'
  #' This function provides ggplot2 theme settings optimized for manuscript figures, featuring no legends
  #' and a bold title using Arial font. It also sets the plot to use void theme which removes most plot background elements.
  #'
  #' @return A ggplot2 theme object tailored for use in manuscript figures.
  #' @import ggplot2
  manuscript_plot_settings <- theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(hjust=0.5, face='bold', family="Arial")) + 
    xlab()
}



## |||||||||||||
## Plot Components ----
## |||||||||||||
# Plots the functional components
plot.components <- function(locations, loadings, type = "points",
                            size = 1,  colormap = "D", limits = NULL,
                            ncol = 5) {
  #' Plot functional components as points or tiles
  #'
  #' Plots a series of functional components either as points or tiles based on the type parameter.
  #' Supports multiple plotting through arranging sub-plots.
  #'
  #' @param locations Data frame with coordinates for plotting.
  #' @param loadings Matrix with component loadings.
  #' @param type Character, either "points" or "tile" indicating the type of plot.
  #' @param size Numeric, size of points or tiles.
  #' @param colormap String, the name of the colormap to use.
  #' @param limits Numeric vector, limits for data scaling in color mapping.
  #' @param ncol Integer, number of columns for the layout of sub-plots.
  #' @return A combined ggplot object with multiple plots of components.
  #' @import ggplot2
  #' @import viridis
  
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
  #' Guide plot for optimal number of components selection
  #'
  #' Generates plots to aid in the selection of an optimal number of components based on normalized residuals.
  #' The plot includes a line and point plot showing the percentage of variance explained and gains with each component.
  #'
  #' @param residuals_norm Numeric vector of normalized residuals.
  #' @param nComp Integer, total number of components considered.
  #' @param nComp_opt Integer, optimal number of components.
  #' @return A plot object with lines and points indicating the performance metrics per number of components.
  #' @import ggplot2
  
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
  
  ### TESTING WITH PIETRO
  print(data_plot)
  ###TESTING WITH PIETRO
  
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


plot.Scomps <- function(locations, loadings,
                        size = 1, colormap_continuous = coolwarm(200),
                        colormap_discrete = "Set2", discrete = FALSE, 
                        ncol = 5, title_prefix = "Comp",
                        extra_formatting = NULL, is_equal_scale = TRUE,
                        custom_colors=NULL, custom_titles=NULL) {
  #' Plot series of component loadings as point plots. Allows scaling normalization based 
  #' on 1st component limits
  #'
  #' Generates multiple point plots for each set of component loadings provided, supporting both
  #' continuous and discrete color scales. Intended for visualization of PCA-like loadings.
  #'
  #' @param locations A data frame or matrix of x and y coordinates.
  #' @param loadings A matrix or data frame where each column represents component loadings.
  #' @param size Numeric, size of the points in the plot.
  #' @param colormap_continuous A color palette for continuous data values.
  #' @param colormap_discrete A color palette name for discrete data values.
  #' @param discrete Logical, indicating if color mapping should be discrete.
  #' @param ncol Integer, number of columns for arranging the plots.
  #' @param title_prefix String, prefix for the titles of individual plots.
  #' @param extra_formatting A function that returns additional ggplot2 theme settings.
  #' @param is_equal_scale Logical, whether to use the same scale across all plots.
  #' @param custom_colors A vector of colors to manually specify for discrete scales.
  #' @param custom_titles A vector of custom titles, one for each plot.
  #' @return A list of ggplot objects or a combined plot object.
  #' @import ggplot2
  
  # Calculate limits based on the first component if is_equal_scale is TRUE.
  loadings <- loadings
  limits <- NULL
  
  plots <- list()
  for(i in 1:ncol(loadings)){
    # If is_equal_scale is FALSE, calculate limits for each component individually.
    if(!is_equal_scale) {
      limits <- range(loadings[, i], na.rm = TRUE)
    } else if (is_equal_scale){
      ###NEED to rescale loadings to make sure it is in the range of the first component
      rangei <- range(loadings[ ,i])
      middle <- (rangei[2] + rangei[1]) /2
      loadings[ ,i] <- loadings[ ,i] - middle
      if (i == 1){
        limits <- range(loadings[, 1], na.rm = TRUE)
      }
    }
    
    #Set title
    title = paste(title_prefix, i)
    if (!is.null(custom_titles)){
      if (length(custom_titles) != ncol(loadings)){
        stop("Error: The argument 'custom_titles' is not the same length of the number of subplots ('loadings') in 
             plot.Scomps")
      }
      title = custom_titles[i]
    }
    
    #Plot
    plots[[i]] <- plot.points(locations, loadings[,i], size = size,
                              limits = limits, cmap_continuous = colormap_continuous,
                              cmap_discrete = colormap_discrete, discrete = discrete,
                              title = title, custom_colors = custom_colors)
    
    # Apply extra formatting if provided
    if (!is.null(extra_formatting)) {
      extra_theme <- extra_formatting()
      plots[[i]] <- plots[[i]] + extra_theme
    }
  }
  
  nComp <- length(plots) # Assuming that nComp should be the number of plots
  plot <- arrangeGrob(grobs = plots, ncol = min(ncol, nComp))
  
  return(plot)
}







plot.Scomps2 <- function(locations, loadings,
                         size = 1, colormap_continuous = coolwarm(200),
                         colormap_discrete = "Set2", discrete = FALSE, 
                         ncol = 5, title_prefix = "Comp",
                         extra_formatting = NULL, is_equal_scale = TRUE,
                         custom_colors=NULL, custom_titles=NULL) {
  #' Plot series of component loadings as point plots. This version has an adjustment for scaling the
  #' components based on a global limit between all components as compared to jsut the 1st component
  #'
  #' Generates multiple point plots for each set of component loadings provided, supporting both
  #' continuous and discrete color scales. Intended for visualization of PCA-like loadings.
  #'
  #' @param locations A data frame or matrix of x and y coordinates.
  #' @param loadings A matrix or data frame where each column represents component loadings.
  #' @param size Numeric, size of the points in the plot.
  #' @param colormap_continuous A color palette for continuous data values.
  #' @param colormap_discrete A color palette name for discrete data values.
  #' @param discrete Logical, indicating if color mapping should be discrete.
  #' @param ncol Integer, number of columns for arranging the plots.
  #' @param title_prefix String, prefix for the titles of individual plots.
  #' @param extra_formatting A function that returns additional ggplot2 theme settings.
  #' @param is_equal_scale Logical, whether to use the same scale across all plots.
  #' @param custom_colors A vector of colors to manually specify for discrete scales.
  #' @param custom_titles A vector of custom titles, one for each plot.
  #' @return A list of ggplot objects or a combined plot object.
  #' @import ggplot2
  
  # Initialize
  loadings <- loadings
  limits <- NULL
  plots <- list()
  
  # Calculate global limits based on all components if is_equal_scale is TRUE
  if(is_equal_scale) {
    all_limits <- range(loadings, na.rm = TRUE)
  }
  
  # Then use 'all_limits' for every plot instead of recalculating per component
  for(i in 1:ncol(loadings)){
    limits <- if(is_equal_scale) all_limits else range(loadings[, i], na.rm = TRUE)
    
    
    #Set title
    title = paste(title_prefix, i)
    if (!is.null(custom_titles)){
      if (length(custom_titles) != ncol(loadings)){
        stop("Error: The argument 'custom_titles' is not the same length of the number of subplots ('loadings') in 
             plot.Scomps")
      }
      title = custom_titles[i]
    }
    
    #Plot
    plots[[i]] <- plot.points(locations, loadings[,i], size = size,
                              limits = limits, cmap_continuous = colormap_continuous,
                              cmap_discrete = colormap_discrete, discrete = discrete,
                              title = title, custom_colors = custom_colors)
    
    # Apply extra formatting if provided
    if (!is.null(extra_formatting)) {
      extra_theme <- extra_formatting()
      plots[[i]] <- plots[[i]] + extra_theme
    }
  }
  
  nComp <- length(plots) # Assuming that nComp should be the number of plots
  plot <- arrangeGrob(grobs = plots, ncol = min(ncol, nComp))
  
  return(plot)
}





## |||||||||||||
## Plot Fields ----
## ||||||||||||||||

# Points version
plot.field_points <- function(locations, data,
                              size = 1, limits = NULL, colormap = "D",
                              discrete = FALSE) {
  #' Plot field data as points on a map
  #'
  #' Plots a field of data points on a coordinate system using ggplot2. Allows coloring points by data values,
  #' with support for both continuous and discrete scales.
  #'
  #' @param locations A data frame or matrix with x and y coordinates.
  #' @param data A vector of data values associated with each location.
  #' @param size Numeric, size of the points.
  #' @param limits Numeric vector of length 2, specifying the range of data values for color mapping.
  #' @param colormap String, name of the colormap from viridis.
  #' @param discrete Logical, whether to use a discrete color scale.
  #' @return A ggplot object representing the data as points with specified aesthetics.
  #' @import ggplot2
  #' @import viridis
  data_plot <- data.frame(locations, value = data)
  colnames(data_plot) <- c("x", "y", "value")
  
  plot <- ggplot() +
    theme_minimal() +
    geom_point(data = data_plot, aes(x = x, y = y, color = value), size = size) +
    coord_fixed() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    standard_plot_settings()
  
  
  if(!discrete){
    if(is.null(limits)){
      plot <- plot + scale_color_viridis(option = colormap)
    } else {
      plot <- plot + scale_color_viridis(option = colormap,
                                         limits = limits)
    }
  } else {
    plot <- plot + scale_color_viridis_d(option = colormap)
  }
  
  return(plot)
}



# Tile version
plot.field_tile <- function(nodes, f,
                            limits = NULL, colormap = "D") {
  #' Plot field data as tiles on a map
  #'
  #' Creates a tile plot of field data using ggplot2, filling tiles based on data values.
  #' Supports both continuous and specified data value ranges for coloring.
  #'
  #' @param nodes A data frame or matrix with x and y coordinates.
  #' @param f A vector of data values to fill tiles.
  #' @param limits Numeric vector of length 2, specifying the range of data values for fill mapping.
  #' @param colormap String, name of the colormap from viridis to apply.
  #' @return A ggplot object representing the data as tiled areas with specified aesthetics.
  #' @import ggplot2
  #' @import viridis
  
  data_plot <- data.frame(nodes, value = f)
  colnames(data_plot) <- c("x", "y", "value")
  
  data_plot <- na.omit(data_plot)
  
  plot <- ggplot() +
    theme_minimal() +
    geom_tile(data = data_plot, aes(x = x, y = y, fill = value)) +
    coord_fixed() +
    standard_plot_settings()
  
  if(is.null(limits)){
    plot <- plot + scale_fill_viridis(option = colormap)
  } else {
    plot <- plot + scale_fill_viridis(option = colormap, 
                                      limits = limits)
  }
  
  return(plot)
}



# Original plot.field_tile
plot.field_tile_original <- function(nodes, f, boundary = NULL,
                            limits = NULL, breaks = NULL, colormap = "D",
                            discrete = FALSE, ISOLINES = FALSE, LEGEND = FALSE) {
  
  if(is.null(f)) {
    return(ggplot() + theme_void())
  }
  
  data <- data.frame(nodes, value = f)
  colnames(data) <- c("x", "y", "value")
  
  data <- na.omit(data)
  
  plot <- ggplot() +
    geom_tile(data = data, aes(x = x, y = y, fill = value)) +
    coord_fixed()
  
  if(!is.null(breaks) || ISOLINES) {
    color = "black"
    limits_real <- range(data$value)
    if(is.null(breaks)) {
      breaks <- seq(limits_real[1], limits_real[2], length = 10)
    }
    breaks_initial <- breaks
    h = breaks[2]-breaks[1]
    if(limits_real[1] < min(breaks)) {
      breaks <- c(sort(seq(min(breaks), limits_real[1]-h, by = -h)[-1]), breaks)
    }
    if(limits_real[2] > max(breaks)) {
      breaks <- c(breaks, seq(max(breaks), limits_real[2]+h, by = h)[-1])
    }
    if(length(breaks) > 2*length(breaks_initial) ) {
      breaks <- breaks_initial
      color = "red"
    }
    plot <- plot +
      geom_contour(data = data, aes(x = x, y = y, z = value), color = color, breaks = breaks)
  }
  
  if (!discrete) {
    if(!is.null(breaks)) {
      h = breaks[2]-breaks[1]
      limits <- limits + c(-h, h)
    }
    if (is.null(limits)) {
      plot <- plot + scale_fill_viridis(
        option = colormap # ,
        # labels = scales::scientific_format()
      )
    } else {
      plot <- plot + scale_fill_viridis(
        option = colormap,
        limits = limits # ,
        # labels = scales::scientific_format()
      )
    }
  } else {
    plot <- plot + scale_fill_viridis_d(
      option = colormap # ,
      # labels = scales::scientific_format()
    )
  }
  
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat), fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(fill = "none")
  }
  
  return(plot)
}




## ||||||||||||||||
## Plot Points ----
## ||||||||||||||||

plot.points <- function(locations, data, size = 1, limits = NULL, 
                        cmap_continuous = coolwarm(200), cmap_discrete="Set2", 
                        discrete = FALSE, title="", custom_colors=NULL) {
  #' Custom point plotting with continuous or discrete color mapping
  #'
  #' This function extends ggplot2 to plot points with options for both continuous and discrete color scales.
  #' It allows for custom color palettes and includes additional formatting options like plot titles.
  #'
  #' @param locations A data frame or matrix of x and y coordinates.
  #' @param data A vector of data values.
  #' @param size Numeric, size of the points.
  #' @param limits Numeric vector of length 2, specifying the range of data values for color mapping.
  #' @param cmap_continuous A color palette for continuous data values.
  #' @param cmap_discrete A color palette name for discrete data values.
  #' @param discrete Logical, whether to use a discrete color scale.
  #' @param title String, title of the plot.
  #' @param custom_colors A vector of colors to manually specify for discrete scales.
  #' @return A ggplot object representing the data as points with specified aesthetics.
  #' @import ggplot2
  
  data_plot <- data.frame(locations, value = data)
  colnames(data_plot) <- c("x", "y", "value")
  
  if(discrete){
    data_plot$value <- as.factor(data_plot$value)
  }
  
  plot <- ggplot() +
    theme_minimal() +
    geom_point(data = data_plot, aes(x = x, y = y, color = value), size = size) +
    coord_fixed() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) 
  
  if(title != ""){
    plot = plot + ggtitle(title)
  }
  
  if(!discrete){
    if(is.null(limits)){
      plot <- plot + scale_color_gradientn(colors = cmap_continuous)
    } else {
      plot <- plot + scale_color_gradientn(colors = cmap_continuous,
                                           limits = limits)
    }
  } else {
    if (is.null(custom_colors)){
      plot <- plot + scale_color_brewer(palette=cmap_discrete)
    } else {
      plot <- plot + scale_color_manual(values=custom_colors)
    }
    
  }
  
  return(plot)
}





## ||||||||||||||||
## Boxplots ----
## ||||||||||||||||
create_boxplot <- function(data, x_var, y_var, fill_var = NULL, fill_colors = NULL,
                           title = NULL, x_label = NULL, y_label = NULL,
                           facet_var = NULL, theme_style = theme_minimal(),
                           save_plot = FALSE, filename = "boxplot.png",
                           plot_width = 8, plot_height = 6, dpi = 300) {
  # Load required packages
  library(ggplot2)
  library(RColorBrewer)  # For color palettes
  
  # Check if variables exist in the dataframe
  vars <- c(x_var, y_var, fill_var, facet_var)
  vars <- vars[!is.null(vars)]  # Remove NULLs
  missing_vars <- vars[!vars %in% names(data)]
  if (length(missing_vars) > 0) {
    stop(paste("Variable(s) not found in the dataframe:", paste(missing_vars, collapse = ", ")))
  }
  
  # Convert variables to factors if they are not already
  data[[x_var]] <- as.factor(data[[x_var]])
  if (!is.null(fill_var)) {
    data[[fill_var]] <- as.factor(data[[fill_var]])
  }
  if (!is.null(facet_var)) {
    data[[facet_var]] <- as.factor(data[[facet_var]])
  }
  
  # Initialize ggplot
  p <- ggplot(data, aes_string(x = x_var, y = y_var))
  
  # Add the boxplot layer with fill if specified
  if (!is.null(fill_var)) {
    p <- p + geom_boxplot(aes_string(fill = fill_var))
  } else {
    p <- p + geom_boxplot()
  }
  
  # Apply custom fill colors if provided
  if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values = fill_colors)
  } else if (!is.null(fill_var)) {
    # Generate a color palette automatically
    num_colors <- length(levels(data[[fill_var]]))
    palette <- brewer.pal(min(num_colors, 8), "Set2")
    if (num_colors > 8) {
      palette <- colorRampPalette(palette)(num_colors)
    }
    p <- p + scale_fill_manual(values = palette)
  }
  
  # Add facets if specified
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste('~', facet_var)))
  }
  
  # Add labels
  p <- p + labs(title = title, x = x_label, y = y_label, fill = fill_var)
  
  # Apply theme
  p <- p + theme_style +
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 12)
    )
  
  # Save the plot if requested
  if (save_plot) {
    ggsave(filename, plot = p, width = plot_width, height = plot_height, dpi = dpi)
  }
  
  # Return the plot object
  return(p)
}






# Plot utilities ----
# |||||||||||||||||||

## Settings ----
## |||||||||||||

standard_plot_settings <- function(){
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
  manuscript_plot_settings <- theme_void() +
    theme(legend.position="none", 
          plot.title=element_text(hjust=0.5, face='bold', family="Arial")) + 
    xlab()
}


## Plot Fields ----
## ||||||||||||||||

# Points version

plot.field_points <- function(locations, data,
                              size = 1, limits = NULL, colormap = "D",
                              discrete = FALSE) {
  
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



## Drew plotting fxns ----
## ||||||||||||||||

plot.points <- function(locations, data,size = 1, limits = NULL, 
                        cmap_continuous = coolwarm(200), cmap_discrete="Set2", 
                        discrete = FALSE, title="", custom_colors=NULL) {
  
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




plot.Scomps <- function(locations, loadings,
                        size = 1, colormap_continuous = coolwarm(200),
                        colormap_discrete = "Set2", discrete = FALSE, 
                        ncol = 5, title_prefix = "Comp",
                        extra_formatting = NULL, is_equal_scale = TRUE,
                        custom_colors=NULL, custom_titles=NULL) {
  
  # Calculate limits based on the first component if is_equal_scale is TRUE.
  loadings <- loadings
  limits <- NULL
  #if(is_equal_scale) {
  #  limits <- range(loadings[, 1], na.rm = TRUE)
  #}
  
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
    
     #range1 <- range(loadings[ ,1])
     #rangei <- range(loadings[ ,i])
     #if (!(rangei[1] > range1[1])){
     # diff <- range1 - rangei
     # loadings[ ,i] <- loadings[ ,i] + diff
     #} else if(!(rangei[2] < range1[2])){
     # diff <- range1 - rangei
     # loadings[ ,i] <- loadings[ ,i] + diff
     #}
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

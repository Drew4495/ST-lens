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
## ||||||||||||||||||||
# Load libraries ----
## ||||||||||||||||||||

library(ggplot2)
library(grid)
library(gridExtra)
library(Metrics)





#=============================================================================#




## ||||||||||||||||||||
# Define Functions ----
## ||||||||||||||||||||

### Calculate residuals
fPCA_residuals <- function(X, scores, loadings) {
  
  residuals <- c()
  for(h in 1:ncol(scores)){
    residuals[h] <- RMSE(X, scores[,1:h] %*% t(loadings[,1:h]))
  }
  names(residuals) <- as.character(1:ncol(scores))
  
  return(residuals)
}



### Guide plot to decide the optimal number of components
plot.nComp_selection_old <- function(residuals_norm, nComp, nComp_opt) {
  
  # Percentage of data explained
  data_plot <- data.frame(id = 0:nComp,
                          name = factor(names(residuals_norm), levels = names(residuals_norm)),
                          residuals_norm = unlist(residuals_norm))
  plot1 <- ggplot(data = data_plot) +
    geom_hline(yintercept = min(data_plot$residuals_norm), linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_line(aes(x = id, y = residuals_norm, group = 1)) +
    geom_point(aes(x = id, y = residuals_norm), size = 2) +
    geom_point(aes(x = id[1+nComp_opt], y = residuals_norm[1+nComp_opt]), size = 2, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Residuals") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings()
  
  # Information gained
  data_plot <- data.frame(id = 0:nComp,
                          name = factor(names(residuals_norm), levels = names(residuals_norm)),
                          gain = - unlist(residuals_norm) + c(NaN, unlist(residuals_norm)[1:(nComp)]))
  plot2 <- ggplot(data = data_plot) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
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












### Guide plot to decide the optimal number of components
plot.nComp_selection <- function(residuals_norm, nComp, nComp_opt, theme_custom = NULL) {
  
  # Percentage of data explained
  data_plot <- data.frame(
    id = 0:nComp,
    name = factor(names(residuals_norm), levels = names(residuals_norm)),
    residuals_norm = unlist(residuals_norm)
  )
  
  plot1 <- ggplot(data = data_plot) +
    geom_hline(yintercept = min(data_plot$residuals_norm), linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_line(aes(x = id, y = residuals_norm, group = 1)) +
    geom_point(aes(x = id, y = residuals_norm), size = 2) +
    geom_point(aes(x = id[1 + nComp_opt], y = residuals_norm[1 + nComp_opt]), size = 2, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Residuals") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings()
  
  # Apply custom theme if provided
  if (!is.null(theme_custom)) {
    plot1 <- plot1 + theme_custom
  }
  
  # Information gained
  data_plot <- data.frame(
    id = 0:nComp,
    name = factor(names(residuals_norm), levels = names(residuals_norm)),
    gain = -unlist(residuals_norm) + c(NaN, unlist(residuals_norm)[1:nComp])
  )
  
  plot2 <- ggplot(data = data_plot) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_line(aes(x = id, y = gain, group = 1)) +
    geom_point(aes(x = id, y = gain), size = 2) +
    geom_point(aes(x = id[1 + nComp_opt], y = gain[1 + nComp_opt]), size = 2, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Percentage gain") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings()
  
  # Apply custom theme if provided
  if (!is.null(theme_custom)) {
    plot2 <- plot2 + theme_custom
  }
  
  # Title
  title <- textGrob(
    "Optimal nComp selection",
    gp = gpar(fontsize = 18, fontface = 'bold')
  )
  
  # Plot
  plot <- arrangeGrob(title, plot1, plot2, nrow = 3, heights = c(1, 4, 4))
  
  return(plot)
}

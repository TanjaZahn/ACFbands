#' @title Helper function: Plot analytical examples
#'
#' @description This helper function is used to create plots for the analytical examples.
#'
#' @param df data frame
#' @param yvar character string specifying the y-variable.
#' @param colorvar character string specifying the group variable.
#' @param ylabel label of the y-axis.
#' @param hline numeric and optional. If specified, a horizontal line is drawn at that number.
#'
#' @return a plot
#' 


# Function to plot examples
plot_example <- function(df, yvar, colorvar, ylabel, hline = NULL){
  
  yvar <- ensym(yvar)
  colorvar <- ensym(colorvar)
  
  p <- ggplot(df, aes(x = H, y = !!yvar)) +
    geom_line(aes(color = !!colorvar)) +
    scale_color_manual(labels = mylabels,  values = mycolors) +
    labs(x = "H", y = ylabel) +
    mytheme
  
  if(!is.null(hline)){p <-  p + geom_hline(aes(yintercept = hline))}
  
  return(p)
  
}
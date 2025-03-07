#' @title Helper function: Plot multiple inference bands
#'
#' @description This helper function is used in the applications to plot several inference bands in one graph using different colors and (possibly) different linetypes.
#'
#' @param df data frame
#' @param segment logical. Determines whether sample autocorrelations are plotted are vertically  (`TRUE`) or horizontally (`FALSE`). The first is used for significance bands, the latter for confidence bands.
#' @param line_vec vector containing the linetypes of the different inference bands.
#' @param color_vec vector containing the colors of the different inference bands.
#'
#' @return a plot
#' 



# Plot several inference bands in one graph using different colors and (possibly) linetypes
make_plot_color <- function(df, segment, line_vec, color_vec){
  
  # lock in factor level order
  df$type <- factor(df$type, levels = unique(df$type))
  
  # Show autocorrelation as vertical bars
  if(segment == TRUE){
    myplot <-  ggplot(df, aes(x = h)) + geom_segment(aes(y = rho_hat, xend = h, yend = 0))}
  
  # Show autocorrelation as line
  if(segment == FALSE){
    myplot <-  ggplot(df, aes(x = h)) + geom_line(mapping = aes(y = rho_hat), color = "black") }
  
  # Add upper and lower bounds
  myplot <- myplot +
    geom_hline(aes(yintercept = 0), color = "cornsilk4") +
    geom_line(aes(y = lb, color = type , linetype = type, color = type)) +
    geom_line(aes(y = ub, color = type, linetype = type, color = type)) +
    labs(x = "Lag") +
    scale_x_continuous(limits = c(0, H), breaks = seq(0, H, 5)) +
    scale_linetype_manual(values = line_vec) + 
    scale_color_manual(values= color_vec) +
    mytheme +
    scale_y_continuous(limits = c(-1, 1), breaks = round(seq(-1, 1, 0.2), 1)) +
    theme(legend.position = c(0.1, 0.14), legend.margin=margin(c(0,1,1,1)))
  
  myplot
  
}
#' @title Plot autocorrelations with inference bands.
#'
#' @description This function is used in the three main functions to create a plot. 
#'
#' @param rho_hat A vector of length `H` containing the sample autocorrelations.
#' @param lb A vector of length `H` containing the lower bound of the inference band.
#' @param ub A vector of length `H` containing the upper bound of the inference band.
#' @param color color of the lines.
#' @param segment logical. Determines whether sample autocorrelations are plotted are vertically  (`TRUE`) or horizontally (`FALSE`). The first is used for significance bands, the latter for confidence bands.

#'
#' @return a vector of length `H` containing the estimated autocorrelations at each lag.
#' 
#' @export


make_plot <- function(rho_hat, lb, ub, color, segment){
  
  df <- data.frame(h = 1:length(rho_hat), 
                   rho_hat = rho_hat,
                   lb =  lb,
                   ub =  ub)
  
  # Show autocorrelation as vertical bars
  if(segment == TRUE){
    myplot <-  ggplot(df, aes(x = h)) + geom_segment(aes(y = rho_hat, xend = h, yend = 0))}
  
  # Show autocorrelation as line
  if(segment == FALSE){
    myplot <-  ggplot(df, aes(x = h)) + geom_line(mapping = aes(y = rho_hat), color = "black")}
  
  # Add upper and lower bounds
  myplot <- myplot +
    geom_hline(aes(yintercept = 0), color = "cornsilk4") +
    geom_line(aes(y = lb), color = color , linetype = "dashed") +
    geom_line(aes(y = ub), color = color, linetype = "dashed") +
    labs(x = "Lag", y = "ACF") +
    theme_bw()
  
}
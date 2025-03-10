#' @title Estimate autocorrelations
#'
#' @description Estimate autocorrelations of a given time series for lags `1,...,H`.
#'
#'
#' @param y a vector of observations
#' @param H maximum number of lags for the ACF
#'
#' @return a vector of length `H` containing the estimated autocorrelations at each lag.
#' 
#' @export



autocor <- function(y, H){
  
  sapply(1:H, function(h) sum((y-mean(y))*(dplyr::lead(y, h) - mean(y)), na.rm = TRUE)/sum((y - mean(y))^2, na.rm = TRUE))
  
}

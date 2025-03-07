#' @title DGP: AR(1)
#'
#' @description Generate data according to an AR(1) process.
#'
#'
#' @param N sample size
#' @param I length of burn-in period. Default is 50.
#' @param phi slope coefficient
#'
#' @return vector of length `N` containing observations generated from an AR(1) process
#' 
#' 


sim_ar1_y <- function(N , I = 50, phi){
  
  # Draw normal error
  epsilon <- rnorm(n = N + I, mean = 0, sd = 1)
  
  # Initialization
  y <- epsilon[1]
  
  # AR model
  for(tt in 2:(N + I)) y[tt] <- phi*y[tt-1] + epsilon[tt]
  
  # Disregard burn-in observations
  y <- y[-(1:I)]
  
}
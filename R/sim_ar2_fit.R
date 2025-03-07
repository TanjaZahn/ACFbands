#' @title Estimate AR(1) model on data from AR(2) DGP
#'
#' @description Generate data according to an AR(2) process and fit an AR(1) model to it.
#'
#'
#' @param N sample size
#' @param I length of burn-in period. Default is 50.
#' @param a intercept
#' @param phi1 coefficient of lag 1
#' @param phi2 coefficient of lag 2
#'
#' @return a fitted model object
#' 
#' 

sim_ar2_fit <- function(N , I = 50, a, phi1, phi2){
  
  # Draw normal error
  epsilon <- rnorm(n = N + I, mean = 0, sd = 1)

  # Initialization
  y <- epsilon[1:2]
  
  # AR(2) model
  for(tt in 3:(N + I)) y[tt] <- a + phi1*y[tt-1] + phi2*y[tt-2] + epsilon[tt]
  
  # Collect results in a data frame
  df <- data.frame(y = y, L1_y = dplyr::lag(y, 1))
  
  # Disregard burn-in observations
  df <- df[-(1:I), ]
  
  # Fit an AR(1) model to the data
  fit <- lm(y ~ L1_y, df)
  
}
#' @title Estimate Bartlett covariance matrix of autocorrelations
#'
#' @description Estimate covariance matrix of autocorrelations according to Bartlett's formula (Bartlett, 1946) using the algorithm proposed by Mélard and Roy (1987).
#'
#'
#' @param y a vector of observations
#' @param H maximum number of lags for the ACF
#' @param L bandwidth used for covariance estimation.
#' @param diag_only logical. If `FALSE` (default), the entire covariance matrix is estimated. Otherwise, only the diagonal elements are estimated.
#' 
#'
#' @return Either the covariance matrix of dimension `H x H`(default) or vector of length `H`containing the diagonal elements of the covariance matrix if `diag_only`has been set to `TRUE`.


covar_bartlett <- function(y, H, L, diag_only = FALSE){
  
  # Define kernel
  w <- function(h){
    
    x <- h/L
    
    if(abs(x) <= 1)
      out <- 1 - abs(x)
    
    if(abs(x) > 1)
      out <- 0
    
    return(out)
  }
  
  # Number of observations
  N <- length(y)
  
  # Estimate autocovariance: 
  gamma <- sapply(0:(N-1), function(h) (1/N)*sum((y-mean(y))*(dplyr::lead(y, h) - mean(y)), na.rm = TRUE))
  
  # Collect autocovariances (two-sided)
  c_k <- c(rev(gamma[-1]), gamma)
  
  # Calculate weights
  w_k <- sapply((-N+1):(N-1), function(i) w(i))
  
  # Weighted sum
  lambda <- sapply(0:(N-1), function(i) sum(w_k*dplyr::lead(w_k, i)*c_k*dplyr::lead(c_k, i), na.rm = TRUE))
  
  # Add zero for index being equal to N
  lambda <- c(lambda, 0)
  
  # Estimate full covariance matrix
  if(diag_only == FALSE){
    
    # Estimate the long-run covariance matrix
    B_hat <- matrix(NA, nrow = H, ncol = H)
    
    for(i in 1:H){
      for(j in 1:H){
        B_hat[i, j] <- (lambda[i+j+1] + lambda[abs(i-j) + 1] - 2*w(i)*gamma[i + 1]*lambda[j+1]/gamma[0 + 1] -
                          2*w(j)*gamma[j + 1]*lambda[i+1]/gamma[0 + 1] + 
                          2*w(i)*w(j)*gamma[i + 1]*gamma[j + 1]*lambda[0+1]/(gamma[0 + 1]^2))/(gamma[0 + 1]^2)
      }
    }
    
  }
  
  # Estimate only diagonal elements
  if(diag_only == TRUE){
    
    # Estimate the diagonal of long-run covariance matrix
    B_hat <- rep(NA, H)
    
    for(i in 1:H){
      B_hat[i] <- (lambda[i+i+1] + lambda[abs(i-i) + 1] - 2*w(i)*gamma[i + 1]*lambda[i+1]/gamma[0 + 1] -
                          2*w(i)*gamma[i + 1]*lambda[i+1]/gamma[0 + 1] + 
                          2*w(i)*w(i)*gamma[i + 1]*gamma[i + 1]*lambda[0+1]/(gamma[0 + 1]^2))/(gamma[0 + 1]^2)
      
    }
    
  }
  

  # Output
  return(B_hat)
  
}
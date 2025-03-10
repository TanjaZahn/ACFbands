#' @title Equicoordinate quantile
#'
#' @description Calculate the `tau` quantile of `M` given a covariance matrix `Sigma`.
#'
#'
#' @param Sigma covariance matrix of dimension `H x H`.
#' @param tau probability.
#'
#' @return numeric `tau` quantile.
#' 
#' @export
#' 

equi_quantile <- function(Sigma, tau){
  
  # Check if Sigma is symmetric or has different numbers of rows and columns
  if(nrow(Sigma) != ncol(Sigma)) stop("Sigma has different numbers of rows and columns.")
  
  # Dimension
  H <- nrow(Sigma)
  
  # Vector of expected values
  mu <- rep(0, H)
  
  # Rescale Sigma
  A <- matrix(0, nrow = H, ncol = H)
  diag(A) <- 1/sqrt(diag(Sigma)) 
  
  # New Covariance matrix
  Sigma_new <- A %*% Sigma %*% t(A)
  
  # Quantile: Sup-t Bands
  q <- qmvnorm(p = tau, tail = "both.tails", mean = mu, sigma = Sigma_new)$quantile
  
  # Output
  return(q)
  
}
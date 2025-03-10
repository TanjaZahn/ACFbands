#' @title Covariance matrix for dynamic regressions
#'
#' @description Estimate the autocovariance matrix of the residuals from a dynamic regression with lagged endogeneous regressors.
#'
#'
#' @param fit a fitted model object
#' @param H maximum number of lags for ACF. The default is `30`.
#' @param error a character string specifying whether conditionally homoskedastic (`"hom"`) or heteroskedastic (`"het"`) errors should be computed.
#'
#' @return A list containing
#' \itemize{
#' \item{`estimate` The estimated `H x H` covariance matrix.}
#' \item{`shrinkage_used` logical. If `TRUE`, some elements of the plug-in estimator have been set equal to the identity matrix to ensure positive definiteness.}
#' \item{`shrinkage_lag` The lag `k`from which on elements of the `H x H` plug-in estimator has been set equal to the identity matrix. If `NA`, no shrinkage was necessary to ensure positive definiteness of the plug-in estimator.}
#' }
#' 
#' @export
#' 


covar_dynamic <- function(fit, H, error){
  
  # Has the algorithm been used?
  shrinkage_used <- FALSE
  shrinkage_lag <- NA
  
  # Extract residuals
  res <- fit$residuals
  
  # Sample size
  N <- length(res)
  
  # Estimate Sample Variance of Residuals
  sigma2_hat <- 1/N*sum(res^2)
  
  # Estimate the autocorrelation of residuals
  rho_hat <- autocor(y = res, H = H)
  
  # Number of regressors
  K <- ncol(fit$model) - 1
  
  # Extract the regressors
  X <- matrix(NA, nrow = N, ncol = K)
  for(kk in 1:K){ X[ , kk] <- fit$model[ , (kk+1)]}
  
  # Calculate mean
  xbar <- colMeans(X)
  
  # Estimate Sigma_x
  aux <- matrix(0, nrow = K, ncol = K) # Auxiliary matrix
  for(tt in 1:N){ aux <- aux + (X[tt, ] - xbar) %*% (t(X[tt, ] - xbar )) } # sum
  Sigma_x <- 1/N*aux
  
  # Generate empty Gamma matrix
  Gamma <- matrix(NA, nrow = H, ncol = K)
  
  # Fill in Gamma matrix
  for(h in (1:H)){
    
    # Alternatively:
    c_h <- rep(NA, K)
    for(kk in 1:K){
      c_h[kk] <- 1/(N)*sum(sapply((h+1):N, function(tt) X[tt, kk]*res[tt - h]))
    }
    
    Gamma[h, ] <- c_h
    
  }
  
  
  
  
  ####### Plug-in estimator: Sigma_rho_tilde (possibly negative definite) ----
  
  
  # Conditional Homoskedasticity
  if(error == "hom"){
    
    Sigma_rho_tilde <- diag(H) - (Gamma %*% solve(Sigma_x) %*%  t(Gamma))/sigma2_hat
    
  }
  
  
  # Conditional Heteroskedasticity
  if(error == "het"){
    
    aux <- matrix(0, nrow = K, ncol = K)
    for(tt in 1:N){ aux <- aux + res[tt]^2*(X[tt, ] - xbar) %*% (t(X[tt, ] - xbar ))}
    Sigma_xeps <- 1/N*aux
    
    Sigma_rho_tilde <- diag(H) - 2*(Gamma %*% solve(Sigma_x) %*%  t(Gamma))/sigma2_hat + 
      (Gamma %*% solve(Sigma_x) %*%  Sigma_xeps %*% solve(Sigma_x) %*% t(Gamma))/((sigma2_hat)^2)
    
  }
  
  
  ##### Algorithm: Positive definite Sigma_rho_hat -----------------------------
  
  # Set estimator equal to identity matrix
  Sigma_rho_hat <- diag(H)
  
  # Check positive definiteness of plug-in estimator and replace elements
  k <- 1
  while(k <= H){
    
    # Calculate the eigenvalues of the matrix
    eigenvalues <- eigen(Sigma_rho_tilde[(1:k), (1:k)])$values
    
    # If eigenvalues are positive
    if(all(eigenvalues > 0)){ 
      
      # Replace entries
      Sigma_rho_hat[1:k, 1:k] <- Sigma_rho_tilde[1:k, 1:k]
      
      k <- k+1 } else{
        
        # Keep track that the algorithm has been used    
        shrinkage_used <- TRUE
        
        # At which lag has the shrinkage been used
        shrinkage_lag <- k
        
        # Stop the algorithm
        k <- H+1
        
      }
  }
  
  # Output ---------------------------------------------------------------------
  
  out <- list(estimate = Sigma_rho_hat,
              shrinkage_used = shrinkage_used,
              shrinkage_lag = shrinkage_lag)
  
  
} 


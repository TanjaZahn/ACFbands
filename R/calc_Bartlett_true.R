#' @title Calculate the true covariance matrix for a known AR(1) process
#'
#' @description Calculate true covariance matrix according to the Bartlett formula for an AR(1) process using the formulas provided by Cavazos-Cadena (1994).
#'
#'
#' @param phi coefficient of the AR(1) process
#' @param H maximum number of lags for the ACF
#'
#' @return a `H x H` covariance matrix
#' 




calc_Bartlett_true <- function(phi, H){
  
  # True autocorrelation
  rho <- sapply(1:H, function(h) phi^h)
  
  # Calculate a
  a <- 1 - phi^2
  
  # Calculate R
  R <- sapply(1:H, function(i){
    sapply(1:H, function(j){
      
      if(i == j){
        rho_h <- 1
      } else{
        rho_h <- rho[abs(i - j)]
      }
      
      rho_h
    })
  }) 
  
  # Helper function: Returns the 
  #         - value of the AR(1)-coefficient if index = 1 and 
  #         - 0 otherwise.
  get_arcoef <- function(index){
    
    if(index == 1){
      out <- phi
    }
    
    if(index != 1){
      out <- 0
    }
    
    out
  }
  
  
  # Calculate C
  C <- matrix(NA, ncol = H, nrow = H)
  for(i in 1:H){
    for(j in 1:H){
      C[i, j] <- identical(i, j) - (get_arcoef(i+j) + get_arcoef(i-j))
    }
  }
  
  # Calculate inverse of C
  C_inverse <- solve(C)
  
  # Calculate covariance matrix B (Bartlett)
  B <- a*C_inverse %*% R %*% t(C_inverse)
  
  # Output
  return(B)
  
}

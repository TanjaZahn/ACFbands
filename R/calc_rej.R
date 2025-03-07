#' @title Helper function: Calculate frequency of rejections
#'
#' @description Calculate frequency of rejections for significance bands.
#'
#'
#' @param filename Character string of the filename containing the simulation results.
#' @param add_params Parameters to cbind
#'
#' @return A data frame (one row) containing the avaerage over simulations.
#' 
#' 

calc_rej <- function(filename, add_params){
  
  if(file.exists(filename)){
    
    # Load estimated inference bands
    load(filename)
    
    # Calculate frequency of rejection in simulated samples
    mat_r <- sapply(1:R_mc, function(r){
      
      rho_hat <- estimates_r[[r]]$rho_hat
      sig_band <- estimates_r[[r]]$sig_band
      
      # Reject if rho_hat is outside the bands at h
      rej_h <- (rho_hat < sig_band[ , "lb"]|  rho_hat > sig_band[ , "ub"])
      
      # Reject if rho_hat is outside the bands for at least one h
      rej <- max(rej_h)
      
      # Calculate the average width
      avg_width <- mean(sig_band[ , "width"])
      
      # Combine results
      c(rej = rej, avg_width = avg_width)
      
    })
    
    # Take the mean over R_mc simulations and add parameter values
    out <- cbind(t(rowMeans(mat_r)), add_params)
    
  } else{
    
    # Take the mean over R_mc simulations and add parameter values
    out <- cbind(rej = NA, avg_width = NA, add_params)
    
  }
  
  return(out)
  
}



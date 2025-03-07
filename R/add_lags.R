#' @title Helper function: Create lagged variables 
#'
#' @description Create lagged variables and add it to the data frame.
#'
#' @param df data frame
#' @param colname a character string specifiyng the name of the variable we want to lag.
#' @param lag_vec vector containing the lags that should be created, i.e. c(1, 4) lags the variable by 1 and by 4.
#'
#' @return the original data frame expandend by columns containing the lagged variable.
#' 

add_lags <- function(df, colname, lag_vec){
  
  # Generate lags
  mat_lags <- sapply(lag_vec, function(i){
    
    x <- df %>% pull(any_of(colname))
    dplyr::lag(x, i)
  })
  
  # Name columns
  colnames(mat_lags) <- paste0("L", lag_vec, "_", colname)
  
  # Add lags to data frame
  tibble(df, as_tibble(mat_lags))
  
}
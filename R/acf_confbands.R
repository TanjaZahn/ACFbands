#' @title Autocorrelation estimation with confidence bands (Time Series, Static Regression)
#'
#' @description Estimate and plot autocorrelations with simultaneous confidence bands of level `1-alpha`. Applicable either to a stationary time series or to the residuals from a static regression.
#'
#'
#' @param y a vector of observations or residuals from a static regression, for which the ACF should be computed.
#' @param H maximum number of lags for ACF. The default is `25`.
#' @param type a character string specifying the type of the confidence band to be computed. Allowed are `"sup-t"` (the default), `"bonferroni"` or `"pointwise"`.
#' @param L bandwidth used for covariance estimation. Default is `L = sqrt(length(y))`.
#' @param alpha significance level. Default is `0.05`.
#' @param plot logical. If `TRUE` (default), a plot of the ACF with confidence bands is provided. If `FALSE`, a list of results is provided.
#' @param B_hat an optional `H x H` estimated covariance matrix of the autocorrelations. If not provided (the default), the function will estimate the covariance matrix.
#' @param B_hat_diag an optional vector of length `H` containing the elements of the estimated covariance matrix of the autocorrelations. If not provided (the default), the function will estimate the the necessary elements of the covariance matrix to construct the inference bands. This parameter should only be used in combination with `type` set to `"bonferroni"` or `"pointwise"`. If used with `sup-t`, the input provided to `B_hat_diag` is ignored and the whole covariance matrix is estimated by the function. 
#'
#'
#' @return If `plot = TRUE`, a plot of the ACF with confidence bands is provided. If `plot = FALSE`, a list with the following elements is returned:
#' \itemize{
#' \item{`rho_hat` a vector of length `H` containing the estimated autocorrelations at each lag.}
#' \item{`conf_band` a matrix of dimension `H x 3` containing the lower bound (`lb`), upper bound (`ub`) and the width of the confidence band (`width`).}
#' \item{`B_hat` or `B_hat_diag` estimated `H x H` covariance matrix of autocorrelations for `type = "sup_t"` or vector of length `H` containing the diagonal elements of covariance matrix for `type %in% c("bonferroni", "pointwise")` .}
#' }
#' @examples
#' ## Create artificial data from an AR(1) process
#'  N <- 100 # length of the time series
#'  I <- 50 # burn-in period
#'  phi <- 0.5 # coefficient of first lag in AR(1) process
#'  epsilon <- rnorm(n = N + I, mean = 0, sd = 1) # draw error (I values for burn-in period)
#'  y <- epsilon[1] # Initialization
#'  for(tt in 2:(N + I)) y[tt] <- phi*y[tt-1] + epsilon[tt] # AR model
#'  y <- y[-(1:I)] # disregard burn-in observations
#'  
#'  ## Plot autocorrelations with sup-t confidence bands (default)
#'  acf_confbands(y)
#'  
#'  ## Compare to pointwise confidence bands
#'  acf_confbands(y, type = "pointwise")
#'  
#'  ## Obtain results as a list
#'  res <- acf_confbands(y, plot = FALSE)
#'  res
#'  
#'  
#'
#' @export
#' 


acf_confbands <- function(y, H = 25, type = "sup-t", L = sqrt(length(y)), alpha = 0.05, plot = TRUE, B_hat = NULL, B_hat_diag = NULL){
  
  # Check inputs
  if(!type %in% c("sup-t", "bonferroni", "pointwise")) stop("Please choose a suitable type.")
  if (L <= 0) stop( "Block length L must be positive.")
  if (H <= 0) stop( "Maximum number of lags H must be positive.")
  if (alpha <= 0) stop( "Significance level alpha must be positive.")
  if (!is.logical(plot)) stop( "plot must be logical.")
  
  # Number of observations
  N <- length(y)
  
  # Estimate the autocorrelation from the original sample
  rho_hat <- autocor(y = y, H = H)
  
  # Estimate covariance matrix ------------------------------------------------
  if(type == "sup-t"){
    if(is.null(B_hat)){B_hat <- covar_bartlett(y = y, H = H, L = L)}
    B_hat_diag <- diag(B_hat)
  } 
  
  if(type %in% c("bonferroni", "pointwise")){
    if(is.null(B_hat_diag)){
      if(is.null(B_hat)){
        B_hat_diag <- covar_bartlett(y = y, H = H, L = L, diag_only = TRUE)
      }
      if(!is.null(B_hat)){
        B_hat_diag <- diag(B_hat)
      }
    }

  }
  
  # Calculate quantile ----------------------------------------------------------
  
  # Sup-t bands
  if(type == "sup-t"){q <- equi_quantile(tau = 1 - alpha, Sigma = B_hat)}
  
  # Bonferroni Bands 
  if(type == "bonferroni"){q <- qnorm(p = 1-alpha/(2*H))}
  
  # Pointwise bands
  if(type == "pointwise"){q <- qnorm(p = (1-alpha/2))}
  
  
  # Construct confidence bands -------------------------------------------------
  
  lb <- rho_hat - sqrt(B_hat_diag/N)*q # lower bound
  ub <- rho_hat + sqrt(B_hat_diag/N)*q # upper bound
  width <- ub - lb # width
  
  conf_band <- cbind(lb, ub, width)
  rownames(conf_band) <- paste0("h", 1:H)
  
  # Output -----------------------------------------------------------
  
  if(plot == TRUE){
    
    out <- make_plot(segment = FALSE, color = "chartreuse3",  rho_hat = rho_hat, lb = lb, ub = ub)
    
  }else{
    if(type == "sup-t"){ 
      out <- list(rho_hat = rho_hat, 
                  conf_band = conf_band,
                  B_hat = B_hat) }
    
    if(type %in% c("bonferroni", "pointwise")){ 
      out <- list(rho_hat = rho_hat, 
                  conf_band = conf_band,
                  B_hat_diag = B_hat_diag)}
  }
  
  return(out)
  
}

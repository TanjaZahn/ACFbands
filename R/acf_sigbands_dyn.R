#' @title Autocorrelation estimation with significance bands for dynamic regressions
#'
#' @description Estimate and plot autocorrelations with simultaneous significance bands for the null hypothesis of white noise, i.e. the absence of serial correlation. Applicable to the residuals from a dynamic regression with lagged endogenous regressors. If the sample autocorrelations leave the significance bands at any lag, the null hypothesis is rejected at significance level `alpha`.
#'
#'
#' @param fit a fitted model object
#' @param H maximum number of lags for ACF. The default is `25`.
#' @param type a character string specifying the type of the significance band to be computed. Allowed are `"simultaneous"` (the default) or `"pointwise"`.
#' @param error a character string specifying whether conditionally homoskedastic (`"hom"`) or heteroskedastic (`"het"`) errors should be computed.
#' @param alpha significance level. Default is `0.05`.
#' @param plot logical. If `TRUE` (default), a plot of the ACF with confidence bands is provided. If `FALSE`, a list of results is provided.
#' @param Sigma_rho_hat optional named list containing the estimated covariance matrix under the name `estimate`. If not provided (the default), the function will estimate the covariance matrix.
#'
#' @return If `plot = TRUE`, a plot of the ACF with significance bands is provided. If `plot = FALSE`, a list with the following elements is returned:
#' \itemize{
#' \item{`rho_hat` a vector of length `H` containing the estimated autocorrelations at each lag.}
#' \item{`sig_band` a matrix of dimension `H x 3` containing the lower bound (`lb`), upper bound (`ub`) and the width of the confidence band (`width`).}
#' \item{`Sigma_rho_hat` A list containing
#'        \itemize{
#'          \item `estimate`estimated `H x H` covariance matrix of autocorrelations.}
#'           \item{`shrinkage_used` logical. If `TRUE`, some elements of the plug-in estimator have been set equal to the identity matrix to ensure positive definiteness.}
#'          \item{`shrinkage_lag` The lag `k`from which on elements of the `H x H` plug-in estimator has been set equal to the identity matrix. If `NA`, no shrinkage was necessary to ensure positive definiteness of the plug-in estimator.}
#'        }
#' }
#' @examples
#' ## Create artificial data from an AR(2) process
#'  N <- 100 # length of the time series
#'  I <- 50 # burn-in period
#'  phi1 <- 0.5 # coefficient of first lag in AR(1) process
#'  phi2 <- 0.125 # coefficient of second lag in AR(1) process
#'  epsilon <- rnorm(n = N + I, mean = 0, sd = 1) # draw error
#'  y <- epsilon[1:2] # Initialization
#'  for(tt in 3:(N + I)) y[tt] <- phi1*y[tt-1] + phi2*y[tt-2] + epsilon[tt]# AR model
#'  df <- data.frame(y = y, L1_y = dplyr::lag(y, 1)) # Collect results in a data frame
#'  df <- df[-(1:I), ] # Disregard burn-in observations
#'  
#'  ## Fit an AR(1) model to the data
#'  fit <- lm(y ~ L1_y, df)
#'  
#'  ## Plot autocorrelations of residuals with simultaneous significance bands 
#'  ## assuming conditionally homoskedastic errors (default)
#'  acf_sigbands_dyn(fit)
#'  
#'  ## Compare to conditionally heteroskedastic errors
#'  acf_sigbands_dyn(fit, error = "het")
#'  
#'  ## Compare to pointwise confidence bands with conditionally homoskedastic errors
#'  acf_sigbands_dyn(fit, type = "pointwise", error = "hom")
#'  
#'  ## Obtain results as a list
#'  res <- acf_sigbands_dyn(fit, plot = FALSE)
#'  res
#'  
#'  
#'
#' @export
#' 


acf_sigbands_dyn <- function(fit, H = 25, type = "simultaneous", error = "hom", alpha = 0.05, plot = TRUE, Sigma_rho_hat = NULL){
  
  # Check inputs
  if(!type %in% c("pointwise", "simultaneous")) stop("Please choose a suitable type.")
  if(!error %in% c("hom", "het")) stop("Please choose a suitable type for the error.")
  if (H <= 0) stop( "Maximum number of lags H must be positive.")
  if (alpha <= 0) stop( "Significance level alpha must be positive.")
  if (!is.logical(plot)) stop( "plot must be logical.")
  
  
  # Extract residuals
  res <- fit$residuals
  
  # Sample size
  N <- length(res)
  
  # Estimate the autocorrelation of residuals
  rho_hat <- autocor(y = res, H = H)
  
  # Estimate Sigma_rho_hat
  if(is.null(Sigma_rho_hat)){ Sigma_rho_hat <- covar_dynamic(fit, H, error)}
  
  # Simultaneous bands
  if(type == "simultaneous"){q <-  equi_quantile(tau = 1 - alpha, Sigma = Sigma_rho_hat$estimate)}
  
  # Pointwise bands
  if(type == "pointwise"){q <- qnorm(1-alpha/2)}
  
  # Construct significance bands
  lb <- - sqrt(diag(Sigma_rho_hat$estimate)/N)*q # lower bound
  ub <-+ sqrt(diag(Sigma_rho_hat$estimate)/N)*q # upper bound
  width <- ub - lb # width
  
  sig_band <- cbind(lb, ub, width)
  rownames(sig_band) <- paste0("h", 1:H)
  
  # Output ---------------------------------------------------------------------
  
  if(plot == TRUE){
    
    out <- make_plot(segment = TRUE, color = "blue4", rho_hat = rho_hat, lb = lb, ub = ub)
    
  }else{
    
    out <- list(rho_hat = rho_hat, 
                sig_band = sig_band,
                Sigma_rho_hat = Sigma_rho_hat)
  }
  
  return(out)
  
}

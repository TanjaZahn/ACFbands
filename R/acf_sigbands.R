#' @title Autocorrelation estimation with significance bands (Time Series, Static Regression)
#'
#'@description Estimate and plot autocorrelations with simultaneous significance bands for the null hypothesis of white noise, i.e. the absence of serial correlation. Applicable either to a stationary time series or to the residuals from a static regression. If the sample autocorrelations leave the significance bands at any lag, the null hypothesis is rejected at significance level `alpha`.
#'
#'
#'
#' @param y a vector of observations or residuals from a static regression, for which the ACF should be computed.
#' @param H maximum number of lags for ACF. The default is `25`.
#' @param type a character string specifying the type of the significance band to be computed. Allowed are `"simultaneous"` (the default) or `"pointwise"`.
#' @param alpha significance level. Default is `0.05`.
#' @param plot logical. If `TRUE` (default), a plot of the ACF with confidence bands is provided. If `FALSE`, a list of results is provided.
#'
#' @return If `plot = TRUE`, a plot of the ACF with significance bands is provided. If `plot = FALSE`, a list with the following elements is returned:
#' \itemize{
#' \item{`rho_hat` a vector of length `H` containing the estimated autocorrelations at each lag.}
#' \item{`sig_band` a matrix of dimension `H x 3` containing the lower bound (`lb`), upper bound (`ub`) and the width of the confidence band (`width`).}
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
#'  ## Plot autocorrelations with simultaneous significance bands (default)
#'  acf_sigbands(y)
#'  
#'  ## Compare to pointwise confidence bands
#'  acf_sigbands(y, type = "pointwise")
#'  
#'  ## Obtain results as a list
#'  res <- acf_sigbands(y, plot = FALSE)
#'  res
#'  
#'  
#'
#' @export
#' 


acf_sigbands <- function(y, H = 25, type = "simultaneous",  alpha = 0.05, plot = TRUE){
  
  # Check inputs
  if(!type %in% c("simultaneous", "pointwise")) stop("Please choose a suitable type.")
  if (H <= 0) stop( "Maximum number of lags H must be positive.")
  if (alpha <= 0) stop( "Significance level alpha must be positive.")
  if (!is.logical(plot)) stop( "plot must be logical.")
  
  
  # Estimate the autocorrelation from the original sample
  rho_hat <- autocor(y = y, H = H)
  
  # Sample size
  N <- length(y)
  
  # Quantile for simultaneous bands
  if(type == "simultaneous"){ q <- sqrt(qchisq(p = (1-alpha)^(1/H), df = 1))}
  
  # Quantile for pointwise bands
  if(type == "pointwise"){ q <- qnorm(1-alpha/2)}
  
  # Construct significance bands
  lb <- rep(- sqrt(1/N)*q, H)# lower bound
  ub <- rep(sqrt(1/N)*q, H) # upper bound
  width <- ub - lb # width
  
  sig_band <- cbind(lb, ub, width)
  rownames(sig_band) <- paste0("h", 1:H)
  
  # Output --------------------------------------------------------------------
  
  if(plot == TRUE){
    
    out <- make_plot(segment = TRUE, color = "blue4", H = H, rho_hat = rho_hat, lb = lb, ub = ub)
 
     }else{
    
    out <- list(rho_hat = rho_hat, 
                sig_band = sig_band)
  }
  

  return(out)
  
}

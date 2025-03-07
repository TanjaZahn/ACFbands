# Load the package
library(devtools)
load_all(".")


# Check function for confidence bands -------------------------------------------

set.seed(50)
N <- 100 # length of the time series
I <- 50 # burn-in period
phi1 <- 0.5 # coefficient of first lag in AR(1) process
H <- 30
epsilon <- rnorm(n = N + I, mean = 0, sd = 1) # draw error (I values for burn-in period)
y <- epsilon[1] # Initialization
for(tt in 2:(N + I)) y[tt] <- phi1*y[tt-1] + epsilon[tt] # AR model
y <- y[-(1:I)] # disregard burn-in observations

acf_confbands(y = y, H = H, type = "sup-t", L = sqrt(length(y)), alpha = 0.1, plot = TRUE)

acf_confbands(y = y, H = 5, type = "sup-t", L = sqrt(length(y)), alpha = 0.1, plot = FALSE)

acf_confbands(y = y, H = 5, type = "bonferroni", L = sqrt(length(y)), alpha = 0.1, plot = FALSE)

acf_confbands(y = y, H = 5, type = "pointwise", L = sqrt(length(y)), alpha = 0.1, plot = FALSE)



# Check function for significance bands -------------------------------------------

set.seed(50)
N <- 100 # length of the time series
I <- 50 # burn-in period
phi1 <- 0.5 # coefficient of first lag in AR(1) process
H <- 30
epsilon <- rnorm(n = N + I, mean = 0, sd = 1) # draw error (I values for burn-in period)
y <- epsilon[1] # Initialization
for(tt in 2:(N + I)) y[tt] <- phi1*y[tt-1] + epsilon[tt] # AR model
y <- y[-(1:I)] # disregard burn-in observations

acf_sigbands(y = y, H = H, type = "simultaneous", alpha = 0.1, plot = TRUE)

acf_sigbands(y = y, H = 5, type = "simultaneous", alpha = 0.1, plot = FALSE)

acf_sigbands(y = y, H = 5, type = "pointwise", alpha = 0.1, plot = FALSE)




# Dynamic regressions ----------------------------------------------------------


set.seed(50)
N <- 100 # length of the time series
I <- 50 # burn-in period
phi1 <- 0.5
phi2 <- 0.3
alpha0 <- 1
alpha1 <- 0.3
a <- 0
eta <- rnorm(n = N + I, mean = 0, sd = 1)
h <- eta[1]^2 # take square to avoid negative number
y <- eta[1:2]
for(tt in 2:(N + I)) h[tt] <- alpha0 + alpha1*h[tt-1]*(eta[tt-1])^2
epsilon <- sqrt(h)*eta # alternatively: for(tt in 1:(N + I)) epsilon[tt] <- sqrt(h[tt])*eta[tt]
for(tt in 3:(N + I)) y[tt] <- a + phi1*y[tt-1] + phi2*y[tt-2] + epsilon[tt]
df <- data.frame(y = y, L1_y = dplyr::lag(y, 1))
df <- df[-(1:I), ]
fit <- lm(y ~ L1_y, df)
H <- 30

set.seed(1)
acf_sigbands_dynreg(fit = fit, H = H, type = "simultaneous",  type_error = "const", alpha = 0.1, plot = TRUE)

set.seed(1)
acf_sigbands_dynreg(fit = fit, H = H, type = "simultaneous",  type_error = "HC", alpha = 0.1, plot = TRUE)

set.seed(1)
acf_sigbands_dynreg(fit = fit, H = 5, type = "simultaneous",  type_error = "const", alpha = 0.1, plot = FALSE)

set.seed(1)
acf_sigbands_dynreg(fit = fit, H = 5, type = "simultaneous",  type_error = "HC", alpha = 0.1, plot = FALSE)

set.seed(1)
acf_sigbands_dynreg(fit = fit, H = 5, type = "pointwise",  type_error = "const", alpha = 0.1, plot = FALSE)

set.seed(1)
acf_sigbands_dynreg(fit = fit, H = 5, type = "pointwise",  type_error = "HC", alpha = 0.1, plot = FALSE)






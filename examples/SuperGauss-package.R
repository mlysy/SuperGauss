# Superfast inference for the timescale parameter
# of the exponential autocorrelation function
exp_acf <- function(lambda) exp(-(1:N-1)/lambda)

# simulate data
lambda0 <- 1
N <- 1000
X <- rnormtz(n = 1, acf = exp.acf(lambda0))

# loglikelihood function
Toep <- Toeplitz$new(N) # allocate memory for a Toeplitz matrix object
loglik <- function(lambda) {
  Toep$set_acf(acf = exp.acf(lambda))
  dSnorm(X = X, acf = Toep, log = TRUE)
}

# maximum likelihood estimation
optimize(f = loglik, interval = c(.2, 5), maximum = TRUE)

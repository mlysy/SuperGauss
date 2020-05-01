# Superfast inference for the timescale parameter
# of the exponential autocorrelation function
exp_acf <- function(lambda) exp(-(1:N-1)/lambda)

# simulate data
lambda0 <- 1
N <- 1000
X <- rnormtz(n = 1, acf = exp_acf(lambda0))

# loglikelihood function
# allocate memory for a NormalToeplitz distribution object
NTz <- NormalToeplitz$new(N)
loglik <- function(lambda) {
  NTz$logdens(z = X, acf = exp_acf(lambda))
  ## dSnorm(X = X, acf = Toep, log = TRUE)
}

# maximum likelihood estimation
optimize(f = loglik, interval = c(.2, 5), maximum = TRUE)

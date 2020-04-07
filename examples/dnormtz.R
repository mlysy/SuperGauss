# simulate data
N <- 10 # length of each time series
n <- 3 # number of time series
theta <- 0.1
lambda <- 2
mu <- theta^2 * rep(1, N)
acf <- exp(-lambda * (1:N - 1))

X <- rnormtz(n, acf = acf) + mu

# evaluate log-density
dnormtz(X, mu, acf, log = TRUE)

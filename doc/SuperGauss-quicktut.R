## ----fbmsim, include = FALSE--------------------------------------------------
require(SuperGauss)

N <- 5000 # number of observations
dT <- 1/60 # time between observations (seconds)
H <- .3 # Hurst parameter

tseq <- (0:N)*dT # times at which to sample fBM
npaths <- 5 # number of fBM paths to generate

# to generate fbm, generate its increments, which are stationary
msd <- fbm_msd(tseq = tseq[-1], H = H)
acf <- msd2acf(msd = msd) # convert msd to acf

# superfast method
system.time({
  dX <- rnormtz(n = npaths, acf = acf, fft = TRUE)
})
# fast method (about 3x as slow)
system.time({
  rnormtz(n = npaths, acf = acf, fft = FALSE)
})
# unstructured variance method (much slower)
system.time({
  matrix(rnorm(N*npaths), npaths, N) %*% chol(toeplitz(acf))
})

## ----fbmsim-------------------------------------------------------------------
require(SuperGauss)

N <- 5000 # number of observations
dT <- 1/60 # time between observations (seconds)
H <- .3 # Hurst parameter

tseq <- (0:N)*dT # times at which to sample fBM
npaths <- 5 # number of fBM paths to generate

# to generate fbm, generate its increments, which are stationary
msd <- fbm_msd(tseq = tseq[-1], H = H)
acf <- msd2acf(msd = msd) # convert msd to acf

# superfast method
system.time({
  dX <- rnormtz(n = npaths, acf = acf, fft = TRUE)
})
# fast method (about 3x as slow)
system.time({
  rnormtz(n = npaths, acf = acf, fft = FALSE)
})
# unstructured variance method (much slower)
system.time({
  matrix(rnorm(N*npaths), npaths, N) %*% chol(toeplitz(acf))
})

## ---- fig.width = 10, fig.height = 5, out.width = "90%"-----------------------
# convert increments to position measurements
Xt <- apply(rbind(0, dX), 2, cumsum)

# plot
clrs <- c("black", "red", "blue", "orange", "green2")
par(mar = c(4.1,4.1,.5,.5))
plot(0, type = "n", xlim = range(tseq), ylim = range(Xt),
     xlab = "Time (s)", ylab = "Position (m)")
for(ii in 1:npaths) {
  lines(tseq, Xt[,ii], col = clrs[ii], lwd = 2)
}

## -----------------------------------------------------------------------------
# allocate and assign in one step
Tz <- Toeplitz$new(acf = acf)
Tz

# allocate memory only
Tz <- Toeplitz$new(N = N)
Tz
Tz$set_acf(acf = acf) # assign later

## -----------------------------------------------------------------------------
all(acf == Tz$get_acf()) # extract acf

# matrix multiplication
z <- rnorm(N)
x1 <- toeplitz(acf) %*% z # regular way
x2 <- Tz$prod(z) # with Toeplitz class
x3 <- Tz %*% z # with Toeplitz class overloading the `%*%` operator
range(x1-x2)
range(x2-x3)

# system of equations
y1 <- solve(toeplitz(acf), z) # regular way
y2 <- Tz$solve(z) # with Toeplitz class
y2 <- solve(Tz, z) # same thing but overloading `solve()`
range(y1-y2)

# log-determinant
ld1 <- determinant(toeplitz(acf))$mod
ld2 <- Tz$log_det() # with Toeplitz class
ld2 <- determinant(Tz) # same thing but overloading `determinant()`
                         # note: no $mod
c(ld1, ld2)

## -----------------------------------------------------------------------------
dX <- diff(Xt[,1]) # obtain the increments of a given path
N <- length(dX)

# autocorrelation of fBM increments
fbm_acf <- function(H) {
  msd <- fbm_msd(1:N*dT, H = H)
  msd2acf(msd)
}

# loglikelihood using generalized Schur algorithm
NTz <- NormalToeplitz$new(N = N) # pre-allocate memory
loglik_GS <- function(H) {
  NTz$logdens(z = dX, acf = fbm_acf(H))
}

# loglikelihood using Durbin-Levinson algorithm
loglik_DL <- function(H) {
  dnormtz(X = dX, acf = fbm_acf(H), method = "ltz", log = TRUE)
}

# superfast method
system.time({
  GS_mle <- optimize(loglik_GS, interval = c(.01, .99), maximum = TRUE)
})
# fast method (about 10x slower)
system.time({
  DL_mle <- optimize(loglik_DL, interval = c(.01, .99), maximum = TRUE)
})
c(GS = GS_mle$max, DL = DL_mle$max)

# standard error calculation
require(numDeriv)
Hmle <- GS_mle$max
Hse <- -hessian(func = loglik_GS, x = Hmle) # observed Fisher Information
Hse <- sqrt(1/Hse[1])
c(mle = Hmle, se = Hse)

## -----------------------------------------------------------------------------
T1 <- Toeplitz$new(N = N)
T2 <- T1 # shallow copy: both of these point to the same memory location

# affects both objects
T1$set_acf(fbm_acf(.5))
T1
T2

fbm_logdet <- function(H) {
  T1$set_acf(acf = fbm_acf(H))
  T1$log_det()
}

# affects both objects
fbm_logdet(H = .3)
T1
T2

## -----------------------------------------------------------------------------
T3 <- T1$clone(deep = TRUE)
T1
T3

# only affect T1
fbm_logdet(H = .7)
T1
T3

## -----------------------------------------------------------------------------
# autocorrelation function
exp_acf <- function(t, lambda, sigma) sigma^2 * exp(-abs(t/lambda))
# gradient, returned as a 2-column matrix
exp_acf_grad <- function(t, lambda, sigma) {
  ea <- exp_acf(t, lambda, 1)
  cbind(abs(t)*(sigma/lambda)^2 * ea, # d_acf/d_lambda
        2*sigma * ea) # d_acf/d_sigma
}
# Hessian, returned as an array of size length(t) x 2 x 2
exp_acf_hess <- function(t, lambda, sigma) {
  ea <- exp_acf(t, lambda, 1)
  sl2 <- sigma/lambda^2
  hess <- array(NA, dim = c(length(t), 2, 2))
  hess[,1,1] <- sl2^2*(t^2 - 2*abs(t)*lambda) * ea # d2_acf/d_lambda^2
  hess[,1,2] <- 2*sl2 * abs(t) * ea # d2_acf/(d_lambda d_sigma)
  hess[,2,1] <- hess[,1,2] # d2_acf/(d_sigma d_lambda)
  hess[,2,2] <- 2 * ea # d2_acf/d_sigma^2
  hess
}

# simulate data
lambda <- runif(1, .5, 2)
sigma <- runif(1, .5, 2)
tseq <- (1:N-1)*dT
acf <- exp_acf(t = tseq, lambda = lambda, sigma = sigma)
Xt <- rnormtz(acf = acf)

NTz <- NormalToeplitz$new(N = N) # storage space

# negative loglikelihood function of theta = (lambda, sigma)
# include attributes for gradient and Hessian
exp_negloglik <- function(theta) {
  lambda <- theta[1]
  sigma <- theta[2]
  # acf, its gradient, and Hessian
  acf <- exp_acf(tseq, lambda, sigma)
  dacf <- exp_acf_grad(tseq, lambda, sigma)
  d2acf <- exp_acf_hess(tseq, lambda, sigma)
  # derivatives of NormalToeplitz up to order 2
  derivs <- NTz$hess(z = Xt,
                     dz = matrix(0, N, 2),
                     d2z = array(0, dim = c(N, 2, 2)),
                     acf = acf,
                     dacf = dacf,
                     d2acf = d2acf,
                     full_out = TRUE)
  # negative loglikelihood with derivatives as attributes
  nll <- -1 * derivs$ldens
  attr(nll, "gradient") <- -1 * derivs$grad
  attr(nll, "hessian") <- -1 * derivs$hess
  nll
}

# optimization
system.time({
  mle_fit <- nlm(f = exp_negloglik, p = c(1,1), hessian = TRUE)
})

# display estimates with standard errors
rbind(true = c(lambda = lambda, sigma = sigma),
      est = mle_fit$estimate,
      se = sqrt(diag(solve(mle_fit$hessian))))



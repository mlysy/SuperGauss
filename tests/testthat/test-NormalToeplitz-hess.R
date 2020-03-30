library(SuperGauss)
library(numDeriv)
source("test-functions.R")

context("NormalToeplitz - Hessian matrix.")

nrep <- 10
test_that("The GSchur algorithm returns the correct density", {
  replicate(n = nrep, expr = {
    test_logdens <- function(mu, alpha, X) {
      f <- test_drift_func(mu, N)
      acf <- test_fbm_acf(alpha, dt, N)
      Tz <- Toeplitz(acf = acf)
      dSnorm(X, f, acf = Tz, log = T)
    }
    test_logdens_grad <- function(mu, alpha, X) {
      grad(func = test_logdens, x = mu, alpha = alpha, X = X)
    }
    N <- round(abs(rnorm(n = 1, mean = 200, sd = 5)))
    p <- 2
    dt <- runif(1,0,1)
    alpha <- runif(1,.2,.9)
    mu <- runif(1, 1, 3)
    f <- test_drift_func(mu, N)
    acf <- test_fbm_acf(alpha, dt, N)
    X <- f + rSnorm(n = 1, acf = acf)
    Nt <- NormalToeplitz$new(N = N)
    dz <- cbind(-test_drift_grad(mu, N), rep(0, N))
    dacf <- cbind(rep(0, N), test_fbm_acf_grad(alpha, dt, N))
    d2z <- array(data = 0, dim = c(N, p, p))
    d2z[,1,1] <- -test_drift_hess(mu, N)
    d2acf <- array(data = 0, dim = c(N, p, p))
    d2acf[,2,2] <- test_fbm_acf_hess(alpha, dt, N)
    h11 <- hessian(func = test_logdens, x = mu, alpha = alpha, X = X)
    h12 <- h21 <- grad(func = test_logdens_grad, x = alpha, mu = mu, X = X)
    h22 <- hessian(func = test_logdens, x = alpha, mu = mu, X = X)
    hmat1 <- matrix(c(h11, h12, h21, h22), 2, 2)
    hmat2 <- Nt$hess(z = X - f, dz = dz, d2z = d2z, acf = acf, dacf = dacf, d2acf = d2acf)
    expect_equal(hmat1, hmat2, tolerance = 1e-4)
  })
})


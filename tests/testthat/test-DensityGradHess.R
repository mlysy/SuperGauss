library(SuperGauss)
source("SuperGauss-test-functions.R")
context("density, score and Hessian matrix")

test_that("density, score and Hessian matrix of stationary Gaussian process", {
  N <- round(abs(rnorm(n = 1, mean = 500, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 4, sd = 2)))
  X <- matrix(rnorm(N * d), N, d)
  dT <- 1/60
  theta <- runif(1, min = 0.1, max = 1)
  alpha <- runif(1, min = 0.1, max = 1.9)
  Theta <- c(theta, alpha)
  
  mu <- fbm.mu.SGtest(theta, N)
  acf1 <- fbm.acf.SGtest(alpha, dT, N)
  acf <- Toeplitz(acf = acf1)
  dmu <- dacf <- matrix(0, N, 2)
  dmu[, 1] <- fbm.mu.grad.SGtest(theta, N)
  dacf[, 2] <- fbm.acf.grad.SGtest(alpha, dT, N)
  d2mu <- d2acf <- array(0, c(N, 2, 2))
  d2mu[, 1, 1] <- fbm.mu.hess.SGtest(theta, N)
  d2acf[, 2, 2] <- fbm.acf.hess.SGtest(alpha, dT, N)
  
  expect_equal({
    acf$setAcf(acf1)
    Snorm.grad(X, mu, acf, dmu, dacf)
  }, numDeriv::grad(func = test.grad.hess.SGtest, x = Theta, X = X, acf = acf, dT = dT), 
  tolerance = 1)
  
  expect_equal({
    acf$setAcf(acf1)
    Snorm.Hess(X, mu, acf, dmu, dacf, d2mu, d2acf)
  }, numDeriv::hessian(func = test.grad.hess.SGtest, x = Theta, X = X, acf = acf, dT = dT), 
  tolerance = 1)
})
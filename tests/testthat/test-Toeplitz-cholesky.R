source("SuperGauss-testfunctions.R")

context("Toeplitz - Cholesky Decomposition.")

## Test the Cholesky decomposition results using Levinson's algorithm.

nrep <- 10
test_that("`X = chol(Tz) %*% Z` is computed correctly.", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    Z <- matrix(rnorm(N * d), N, d)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Tz <- toeplitz(acf)
      X1 <- t(chol(Tz)) %*% Z
      X2 <- cholZX(Z = Z, acf = acf)
      expect_equal(X1/max(abs(X1)), X2/max(abs(X2)), tolerance = 1e-6)
    }
  })
})

test_that("`Z = chol(Tz)^{-1} X` is computed correctly.", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    X <- matrix(rnorm(N * d), N, d)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Tz <- toeplitz(acf)
      Z1 <- solve(t(chol(Tz)), X)
      Z2 <- cholXZ(X = X, acf = acf)
      expect_equal(Z1/max(abs(Z1)), Z2/max(abs(Z2)), tolerance = 1e-6)
    }
  })
})

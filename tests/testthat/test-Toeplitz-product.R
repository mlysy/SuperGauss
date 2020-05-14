source("SuperGauss-testfunctions.R")

context("Toeplitz - Product.")

nrep <- 10
test_that("Toeplitz-Matrix multiplication", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    Toep <- Toeplitz$new(N)
    case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern", "zero", "rnd"))
    ncase <- nrow(case.par)
    X <- matrix(rnorm(N * d), N, d)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Toep$set_acf(acf)
      acf.mat <- toeplitz(acf)
      y1 <- Toep %*% X
      y2 <- acf.mat %*% X
      expect_equal(y1, y2, tolerance = 1e-6)
      expect_equal(Toep$prod(X), y2, tolerance = 1e-6)
    }
  })
})

test_that("Matrix-Toeplitz multiplication", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    Toep <- Toeplitz$new(N)
    case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern", "zero", "rnd"))
    ncase <- nrow(case.par)
    X <- matrix(rnorm(N * d), d, N)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Toep$set_acf(acf)
      acf.mat <- toeplitz(acf)
      y1 <- X %*% Toep
      y2 <- X %*% acf.mat
      expect_equal(y1, y2, tolerance = 1e-6)
    }
  })
})

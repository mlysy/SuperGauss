library(SuperGauss)
source("SuperGauss-test-functions.R")
context("Multiply")

test_that("Toeplitz-Matrix multiplication", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- Toeplitz(N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern", "zero", "rnd"),
                          dT = c(1/60, 1/30, 1/15))
  ncase <- nrow(case.par)
  X <- matrix(rnorm(N * d), N, d)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    Toep$setAcf(acf)
    acf.mat <- toeplitz(acf)
    y1 <- Toep %*% X
    y2 <- acf.mat %*% X
    expect_equal(y1, y2, tolerance = 1e-6)
  }
})

test_that("Matrix-Toeplitz multiplication", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- Toeplitz(N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern", "zero", "rnd"),
                          dT = c(1/60, 1/30, 1/15))
  ncase <- nrow(case.par)
  X <- matrix(rnorm(N * d), d, N)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    Toep$setAcf(acf)
    acf.mat <- toeplitz(acf)
    y1 <- X %*% Toep
    y2 <- X %*% acf.mat
    expect_equal(y1, y2, tolerance = 1e-6)
  }
})

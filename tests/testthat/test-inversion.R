library(SuperGauss)
source("SuperGauss-test-functions.R")
context("Inversion")

test_that("Toeplitz-matrix inversion", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- Toeplitz(N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern"),
                          dT = c(1/60, 1/30, 1/15))
  ncase <- nrow(case.par)
  X <- matrix(rnorm(N * d), N, d)
  zero <- matrix(0, N, d)
  Y <- matrix(rnorm(N * d), N, d) * 1e4
  Z <- matrix(rnorm(N * d), N, d) * 1e-4
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    acf.mat <- toeplitz(acf)
    Toep.acf(Toep, acf)
    if(min(eigen(acf.mat)$values) > 0){
      expect_equal(acf.mat %*% Toep.solve(Toep, X), X)
      expect_equal(acf.mat %*% Toep.solve(Toep, zero), zero)
      expect_equal(acf.mat %*% Toep.solve(Toep, Y), Y)
      expect_equal(acf.mat %*% Toep.solve(Toep, Z), Z)
    }
  }
})

library(SuperGauss)
source("SuperGauss-test.functions.R")
context("Multiplication")

test_that("Toeplitz-matrix multiplication", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- new(Toeplitz, N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern", "zero", "rnd"), 
                          dT = c(1/60, 1/30, 1/15), 
                          incr = c(TRUE, FALSE))
  ncase <- nrow(case.par)
  X <- matrix(rnorm(N * d), N, d)
  zero <- matrix(0, N, d)
  Y <- matrix(rnorm(N * d), N, d) * 1e4
  Z <- matrix(rnorm(N * d), N, d) * 1e-4
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    incr <- cp$incr
    acf <- acf.get.SGtest(N, type, dT, incr)
    Toep$AcfInput(acf)
    acf.mat <- toeplitz(acf)
    expect_equal(Toep$Mult(X), acf.mat %*% X)
    expect_equal(Toep$Mult(zero), acf.mat %*% zero)
    expect_equal(Toep$Mult(Y), acf.mat %*% Y)
    expect_equal(Toep$Mult(Z), acf.mat %*% Z)
  }
})
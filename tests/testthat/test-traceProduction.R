library(SuperGauss)
source("SuperGauss-test.functions.R")
context("Multiplication")

test_that("trace of inversion of Toeplitz times Toeplitz", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- new(Toeplitz, N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern"),
                          dT = c(1/60, 1/30, 1/15), 
                          incr = c(TRUE, FALSE))
  ncase <- nrow(case.par)
  acf0 <- rep(0, N)
  acf1 <- rnorm(N)
  acf2 <- rnorm(N) * 1e10
  acf3 <- rnorm(N) * 1e-7
  acf.mat0 <- toeplitz(acf0)
  acf.mat1 <- toeplitz(acf1)
  acf.mat2 <- toeplitz(acf2)
  acf.mat3 <- toeplitz(acf3)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    incr <- cp$incr
    acf <- acf.get.SGtest(N, type, dT, incr)
    acf.mat <- toeplitz(acf)
    Toep$AcfInput(acf)
    expect_equal(Toep$TraceProd(acf0), tr(solve(acf.mat, acf.mat0)))
    expect_equal(Toep$TraceProd(acf1), tr(solve(acf.mat, acf.mat1)))
    expect_equal(Toep$TraceProd(acf2), tr(solve(acf.mat, acf.mat2)))
    expect_equal(Toep$TraceProd(acf3), tr(solve(acf.mat, acf.mat3)))
  }
})
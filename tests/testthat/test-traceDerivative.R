library(SuperGauss)
source("SuperGauss-test-functions.R")
context("derivative of traceProduction")

test_that("derivative of trace of inversion of Toeplitz times Toeplitz", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- Toeplitz(N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern"),
                          dT = c(1/60, 1/30, 1/15))
  ncase <- nrow(case.par)
  acf1 <- rnorm(N)
  acf2 <- rnorm(N)
  acf3 <- rnorm(N)
  acf.mat1 <- toeplitz(acf1)
  acf.mat2 <- toeplitz(acf2)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    acf.mat <- toeplitz(acf)
    Toep$setAcf(acf)
    acf.inv <- solve(acf.mat)
    trace.rst <- trace.SGtest(acf.inv %*% acf.mat1 %*% acf.inv %*% acf.mat2)
    expect_equal(Toep$traceT4(acf1, acf2), trace.rst, tolerance = abs(1e-6 * trace.rst))
  }
})
library(SuperGauss)
source("SuperGauss-test-functions.R")
context("traceProduction")

test_that("trace of inversion of Toeplitz times Toeplitz", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- Toeplitz(N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern"),
                          dT = c(1/60, 1/30, 1/15))
  ncase <- nrow(case.par)
  acf0 <- rep(0, N)
  acf1 <- rnorm(N)
  acf2 <- rnorm(N) * 1e5
  acf3 <- rnorm(N) * 1e-5
  acf.mat0 <- toeplitz(acf0)
  acf.mat1 <- toeplitz(acf1)
  acf.mat2 <- toeplitz(acf2)
  acf.mat3 <- toeplitz(acf3)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    acf.mat <- toeplitz(acf)
    Toep.acf(Toep, acf)
    if(min(eigen(acf.mat)$values) > 0){
      acf.inv <- solve.default(acf.mat)
      expect_equal(traceT2(Toep, acf0), trace(acf.inv %*% acf.mat0))
      expect_equal(traceT2(Toep, acf1), trace(acf.inv %*% acf.mat1))
      expect_equal(traceT2(Toep, acf2), trace(acf.inv %*% acf.mat2))
      expect_equal(traceT2(Toep, acf3), trace(acf.inv %*% acf.mat3))
    }
  }
})
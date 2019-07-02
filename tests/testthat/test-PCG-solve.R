library(SuperGauss)
source("SuperGauss-test-functions.R")
context("PCG-Solve")

test_that("PCG-method inversion", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  P1 <- PCG(N)
  case.par <- expand.grid(type = c("fbm"),
                          dT = c(1/60, 1/30, 1/15), b = c(TRUE, FALSE))
  ncase <- nrow(case.par)
  X <- rnorm(N)
  ntol <- 1e-15
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    acf.mat <- toeplitz(acf)
    if(cp$b) {
      expect_equal(c(acf.mat %*% P1$solve(acf, X, ntol)), X, tolerance = 1e-5)
    } else {
      expect_equal(P1$solve(acf, X, ntol), solve(acf.mat, X), tolerance = 1e-5)
    }
  }
})

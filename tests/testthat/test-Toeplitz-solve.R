library(SuperGauss)
source("test-functions.R")

context("Toeplitz - Solve.")

nrep <- 10
test_that("Toeplitz-matrix inversion", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    Toep <- Toeplitz(N)
    case.par <- expand.grid(type = c("fbm", "matern"), b = c(TRUE, FALSE))
    ncase <- nrow(case.par)
    X <- matrix(rnorm(N * d), N, d)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp$type)
      acf <- test_acf_func(N, type)
      Tmat <- toeplitz(acf)
      Toep$set_acf(acf)
      if(cp$b) {
        expect_equal(Tmat %*% solve(Toep, X), X, tolerance = 1e-5)
      } else {
        expect_equal(Tmat %*% solve(Toep), diag(N), tolerance = 1e-5)
      }
    }
  })
})

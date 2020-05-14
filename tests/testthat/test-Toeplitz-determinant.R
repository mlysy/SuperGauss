source("SuperGauss-testfunctions.R")

context("Toeplitz - Determinant.")

## Test the determinant using Ammar-Gragg's GSchur algorithm.

nrep <- 10
test_that("GSchur algorithm gives correct log-determinant.", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 50, sd = 10)))
    Tz <- Toeplitz$new(N)
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Tmat <- toeplitz(acf)
      Tz$set_acf(acf)
      ldet_R <- determinant(Tmat, logarithm = TRUE)$mod[1]
      expect_equal(determinant(Tz, logarithm = TRUE),
                   ldet_R, tolerance = 1e-6)
      expect_equal(Tz$log_det(), ldet_R, toleranace = 1e-6)
    }
  })
})

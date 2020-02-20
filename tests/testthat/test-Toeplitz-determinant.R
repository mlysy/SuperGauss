library(SuperGauss)
source("test-functions.R")

context("Toeplitz - Determinant.")

## Test the determinant using Ammar-Gragg's GSchur algorithm.

nrep <- 10
test_that("Toeplitz determinant", {
  replicate(n = nrep, expr = {  
    N <- round(abs(rnorm(n = 1, mean = 50, sd = 10)))
    Tz <- Toeplitz(N)
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Tmat <- toeplitz(acf)
      Tz$setAcf(acf)
      expect_equal(determinant(Tz, logarithm = TRUE),
                   determinant(Tmat, logarithm = TRUE)$mod[1],
                   tolerance = 1e-6)
    }
  })
})

library(SuperGauss)
source("test-functions.R")

context("Density")

nrep <- 10
test_that("The GSchur algorithm returns the correct density", {
  replicate(n = nrep, expr = {  
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      X <- rSnorm(n = 1, acf = acf, fft = FALSE)
      Nt <- NormalToeplitz(n = N, p = 1)
      ld1 <- dSnorm(X = X, mu = 0, acf = acf, log = TRUE)
      ld2 <- Nt$logdens(X, acf)
      expect_equal(ld1, ld2)
    }
  })
})


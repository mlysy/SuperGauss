library(SuperGauss)
library(numDeriv)
source("test-functions.R")

context("NormalToeplitz - Full gradient for Auto-Diff algorithms.")

nrep <- 10
test_that("The GSchur algorithm returns the correct density", {
  replicate(n = nrep, expr = {
    test_logdens <- function(X, acf) {
      dSnorm(X, 0, acf = acf, log = T)
    }

    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    p <- 2
    dt <- runif(1,0,1)
    alpha <- runif(1,.2,.9)
    acf <- test_fbm_acf(alpha, dt, N)
    X <- rSnorm(n = 1, acf = acf)
    Nt <- NormalToeplitz(n = N, p = p)
    ans <- Nt$grad_full(X, acf)

    expect_equal(c(jacobian(test_logdens, x = X, acf = acf)), ans$dldz, tolerance = 1e-7)
    expect_equal(c(jacobian(test_logdens, x = acf, X = X)), ans$dldacf, tolerance = 1e-7)
  })
})


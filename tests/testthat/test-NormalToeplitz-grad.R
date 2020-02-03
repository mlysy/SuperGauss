library(SuperGauss)
library(numDeriv)
source("test-functions.R")

context("Density")
nrep <- 10

test_that("The GSchur algorithm returns the correct density", {
  replicate(n = nrep, expr = {  
    test_logdens <- function(mu, alpha, X) {
      f <- test_drift_func(mu, N)
      acf <- test_fbm_acf(alpha, dt, N)
      Tz <- Toeplitz(acf = acf)
      dSnorm(X, f, acf = Tz, log = T)
    }
    N <- round(abs(rnorm(n = 1, mean = 200, sd = 5)))
    p <- 2
    dt <- runif(1, 0, 1)
    mu <- 2
    alpha <- runif(1, .2, .9)
    f <- test_drift_func(mu, N)
    acf <- test_fbm_acf(alpha, dt, N)
    X <- f + rSnorm(n = 1, acf = acf)
    Nt <- NormalToeplitz(n = N, p = 2)
    g1 <- c(grad(func = test_logdens, x = mu, alpha = alpha, X = X), 
            grad(func = test_logdens, x = alpha, mu = mu, X = X))
    dz <- cbind(-test_drift_grad(mu, N), rep(0, N))
    dacf <- cbind(rep(0, N), test_fbm_acf_grad(alpha, dt, N))
    g2 <- Nt$grad(z = X - f, dz = dz, acf = acf, dacf = dacf)

    expect_equal(g1, g2)
  })
})


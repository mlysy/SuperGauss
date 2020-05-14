## library(SuperGauss)
## library(numDeriv)
source("SuperGauss-testfunctions.R")

context("NormalToeplitz - Loglikelihood Gradient.")

nrep <- 10
test_that("NormalToeplitz$grad gives correct result.", {
  replicate(n = nrep, expr = {
    test_logdens <- function(theta, X) {
      mu <- theta[1]
      alpha <- theta[2]
      f <- test_drift_func(mu, N)
      acf <- test_fbm_acf(alpha, dt, N)
      toep_ldens(X - f, acf)
    }
    N <- sample(10:30, 1)
    p <- 2
    dt <- runif(1, 0, 1)
    mu <- 2
    alpha <- runif(1, .2, .9)
    f <- test_drift_func(mu, N)
    acf <- test_fbm_acf(alpha, dt, N)
    X <- f + rnormtz(n = 1, acf = acf)
    Nt <- NormalToeplitz$new(N = N)
    ## g1 <- c(grad(func = test_logdens, x = mu, alpha = alpha, X = X),
    ##         grad(func = test_logdens, x = alpha, mu = mu, X = X))
    g1 <- numDeriv::grad(test_logdens, x = c(mu, alpha), X = X)
    dz <- cbind(-test_drift_grad(mu, N), rep(0, N))
    dacf <- cbind(rep(0, N), test_fbm_acf_grad(alpha, dt, N))
    g2 <- Nt$grad(z = X - f, dz = dz, acf = acf, dacf = dacf)
    expect_equal(g1, g2)
  })
})

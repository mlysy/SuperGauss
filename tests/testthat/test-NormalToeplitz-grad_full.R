library(SuperGauss)
library(numDeriv)
source("test-functions.R")

context("NormalToeplitz - Full gradient for Auto-Diff algorithms.")

nrep <- 10
test_that("The GSchur algorithm returns the correct density", {
  replicate(n = nrep, expr = {
    test_logdens <- function(X, acf) {
      mvtnorm::dmvnorm(X, sigma = toeplitz(acf), log = TRUE)
    }
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    p <- 2
    dt <- runif(1,0,1)
    alpha <- runif(1,.2,.9)
    acf <- test_fbm_acf(alpha, dt, N)
    X <- rnormtz(n = 1, acf = acf)
    Nt <- NormalToeplitz$new(N = N)
    ans_z <- Nt$grad_full(X, acf, calc_dldz = TRUE, calc_dlda = FALSE)
    ans_a <- Nt$grad_full(X, acf, calc_dldz = FALSE, calc_dlda = TRUE)
    ans_za <- Nt$grad_full(X, acf)
    jac_z <- c(jacobian(test_logdens, x = X, acf = acf))
    jac_a <- c(jacobian(test_logdens, x = acf, X = X))
    expect_equal(jac_z, ans_z$dldz, tolerance = 1e-7)
    expect_equal(jac_z, ans_za$dldz, tolerance = 1e-7)
    expect_equal(jac_a, ans_a$dlda, tolerance = 1e-7)
    expect_equal(jac_a, ans_za$dlda, tolerance = 1e-7)
  })
})


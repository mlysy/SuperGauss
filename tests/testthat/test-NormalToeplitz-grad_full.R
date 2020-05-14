## library(SuperGauss)
## library(numDeriv)
source("SuperGauss-testfunctions.R")

context("NormalToeplitz - Loglikelihood Gradient (autodiff).")

nrep <- 10
test_that("NormalToeplitz$grad_full gives correct result.", {
  replicate(n = nrep, expr = {
    ## test_logdens <- function(X, acf) {
    ##   mvtnorm::dmvnorm(X, sigma = toeplitz(acf), log = TRUE)
    ## }
    ## N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    N <- sample(2:10, 1)
    ## p <- 2
    dt <- runif(1,0,1)
    alpha <- runif(1,.2,.9)
    acf <- test_fbm_acf(alpha, dt, N)
    X <- rnormtz(n = 1, acf = acf, fft = FALSE)
    Nt <- NormalToeplitz$new(N = N)
    g1_z <- Nt$grad_full(X, acf, calc_dldz = TRUE, calc_dlda = FALSE)
    g1_a <- Nt$grad_full(X, acf, calc_dldz = FALSE, calc_dlda = TRUE)
    g1_za <- Nt$grad_full(X, acf)
    g2_z <- numDeriv::grad(toep_ldens, x = X, gamma = acf)
    g2_a <- numDeriv::grad(toep_ldens, x = acf, z = X)
    expect_equal(g2_z, g1_z$dldz, tolerance = 1e-7)
    expect_equal(g2_z, g1_za$dldz, tolerance = 1e-7)
    expect_equal(g2_a, g1_a$dlda, tolerance = 1e-7)
    expect_equal(g2_a, g1_za$dlda, tolerance = 1e-7)
  })
})


library(SuperGauss)
source("SuperGauss-testfunctions.R")
context("Toeplitz - Hessian Trace.")

nrep <- 10
test_that("Toeplitz$trace_hess is same in R and C++", {
  replicate(n = nrep, expr = {
    N <- round(runif(1, 20, 40))
    type1 <- sample(c("exp", "fbm", "matern"), 1)
    type2 <- sample(c("exp", "exp2", "fbm", "matern", "zero", "rnd"), 1)
    type3 <- sample(c("exp", "exp2", "fbm", "matern", "zero", "rnd"), 1)
    first20 <- sample(c(TRUE, FALSE), 1)
    first30 <- sample(c(TRUE, FALSE), 1)
    Toep <- Toeplitz$new(N)
    acf <- test_acf_func(N, type1)
    acf2 <- test_acf_func(N, type2, first20)
    acf3 <- test_acf_func(N, type3, first30)
    Toep$set_acf(acf)
    th_R <- solve(toeplitz(acf), toeplitz(acf2))
    th_R <- th_R %*% solve(toeplitz(acf), toeplitz(acf3))
    th_R <- sum(diag(th_R))
    expect_equal(Toep$trace_hess(acf2, acf3), th_R)
  })
})

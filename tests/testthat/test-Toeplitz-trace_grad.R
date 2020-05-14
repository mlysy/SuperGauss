library(SuperGauss)
source("SuperGauss-testfunctions.R")
context("Toeplitz - Gradient Trace.")

nrep <- 10
test_that("Toeplitz$trace_grad is same in R and C++", {
  replicate(n = nrep, expr = {
    N <- round(runif(1, 20, 40))
    type1 <- sample(c("exp", "fbm", "matern"), 1)
    type2 <- sample(c("exp", "exp2", "fbm", "matern", "zero", "rnd"), 1)
    first20 <- sample(c(TRUE, FALSE), 1)
    Toep <- Toeplitz$new(N)
    acf <- test_acf_func(N, type1)
    acf2 <- test_acf_func(N, type2, first20)
    Toep$set_acf(acf)
    tg_R <- solve(toeplitz(acf), toeplitz(acf2))
    tg_R <- sum(diag(tg_R))
    expect_equal(Toep$trace_grad(acf2), tg_R)
  })
})

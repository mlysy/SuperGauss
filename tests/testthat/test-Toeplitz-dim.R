library(SuperGauss)
source("test-functions.R")

context("Toeplitz - Dimensions.")

ntest <- 10
test_that("dim/nrow/ncol", {
  for(ii in 1:ntest){
    N <- sample(20:80, 1)
    acf <- test_acf_func(N, "fbm")
    Tz <- Toeplitz$new(N = N)
    expect_equal(nrow(Tz), N)
    expect_equal(ncol(Tz), N)
    expect_equal(dim(Tz), c(N,N))
  }
})

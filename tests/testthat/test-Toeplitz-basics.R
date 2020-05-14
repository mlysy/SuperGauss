source("SuperGauss-testfunctions.R")

context("Toeplitz - Basics.")

ntest <- 10
test_that("dim/nrow/ncol work as expected.", {
  for(ii in 1:ntest){
    N <- sample(20:80, 1)
    acf <- test_acf_func(N, "fbm")
    Tz <- Toeplitz$new(N = N)
    expect_equal(nrow(Tz), N)
    expect_equal(ncol(Tz), N)
    expect_equal(dim(Tz), c(N,N))
  }
})

ntest <- 5
test_that("shallow and deep clone works as expected", {
  for(ii in 1:ntest){
    N <- sample(5:10, 1)
    acf1 <- exp(-(1:N/N) * runif(1))
    acf2 <- exp(-(1:N/N) * runif(1))
    Tz1 <- Toeplitz$new(acf = acf1)
    # shallow clone
    Tz2 <- Tz1$clone()
    expect_equal(Tz1$get_acf(), acf1)
    expect_equal(Tz2$get_acf(), acf1)
    Tz2$set_acf(acf2)
    expect_equal(Tz1$get_acf(), acf2)
    expect_equal(Tz2$get_acf(), acf2)
    # deep clone
    Tz1$set_acf(acf1)
    Tz2 <- Tz1$clone(deep = TRUE)
    expect_equal(Tz1$get_acf(), acf1)
    expect_equal(Tz2$get_acf(), acf1)
    Tz2$set_acf(acf2)
    expect_equal(Tz1$get_acf(), acf1)
    expect_equal(Tz2$get_acf(), acf2)
  }
})

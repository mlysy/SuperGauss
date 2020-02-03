library(SuperGauss)
source("test-functions.R")

context("Dimensions")

test_that("dim/nrow/ncol", {
  case.par <- expand.grid(N = round(abs(rnorm(n = 10, mean = 50, sd = 10))),
                          pre = c(TRUE, FALSE))
  ncase <- nrow(case.par)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    N <- cp$N
    pre <- cp$pre
    acf <- test_acf_func(N, "fbm")
    if(pre) {
      Tz <- Toeplitz(n = N)
    } else {
      Tz <- Toeplitz(acf = acf)
    }
    expect_equal(nrow(Tz), N)
    expect_equal(ncol(Tz), N)
    expect_equal(dim(Tz), c(N,N))
  }
})

library(SuperGauss)
source("SuperGauss-test.functions.R")
context("Multiplication")

test_that("Toeplitz determinant", {
  N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
  d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
  Toep <- new(Toeplitz, N)
  case.par <- expand.grid(type = c("exp", "exp2", "fbm", "matern"),
                          dT = c(1/60, 1/30, 1/15), 
                          incr = c(TRUE, FALSE))
  ncase <- nrow(case.par)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    incr <- cp$incr
    acf <- acf.get.SGtest(N, type, dT, incr)
    Toep$AcfInput(acf)
    acf.mat <- toeplitz(acf)
    expect_equal(Toep$Det(), log(det(acf.mat)))
  }
})
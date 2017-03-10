library(SuperGauss)
source("SuperGauss-test-functions.R")
context("Determinant")

test_that("Toeplitz determinant", {
  N <- round(abs(rnorm(n = 1, mean = 50, sd = 10)))
  acf <- 1:N
  Toep <- Toeplitz(N)
  case.par <- expand.grid(type = c("exp", "exp2", "matern", "fbm"),
                          dT = c(1/60, 1/30, 1/15))
  ncase <- nrow(case.par)
  for(ii in 1:ncase){
    cp <- case.par[ii, ]
    type <- as.character(cp$type)
    dT <- cp$dT
    acf <- acf.get.SGtest(N, type, dT)
    acf.mat <- toeplitz(acf)
    Toep.acf(Toep, acf)
    if(min(eigen(acf.mat)$values) > 0){
      expect_equal(Toep.det(Toep), log(det(acf.mat)))
    }   
  }
})
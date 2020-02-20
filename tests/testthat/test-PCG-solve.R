library(SuperGauss)
source("test-functions.R")

context("Solve Toeplitz systems using PCG algorithm.")

nrep <- 10
test_that("PCG-method inversion", {
  replicate(n = nrep, expr = {  
    N <- round(abs(rnorm(n = 1, mean = 10, sd = 2)))
    P1 <- PCG(N)
    case.par <- expand.grid(type = c("fbm", "matern"), b = c(TRUE, FALSE))
    ncase <- nrow(case.par)
    X <- rnorm(N)
    ntol <- 1e-15
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp$type)
      acf <- test_acf_func(N, type)
      Tmat <- toeplitz(acf)
      Z <- solve(Tmat, X)
      if(cp$b) {
        expect_equal(max(abs((Tmat %*% P1$solve(acf, X, ntol) - X) / Z)), 0, tolerance = 1e-6)
      } else {
        expect_equal(max(abs((P1$solve(acf, X, ntol)- Z) / Z)), 0, tolerance = 1e-6)
      }
    }
  })  
})

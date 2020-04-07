library(SuperGauss)
source("test-functions.R")

context("Solve Toeplitz systems using PCG algorithm.")

nrep <- 10
test_that("PCG-method inversion", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 10, sd = 2)))
    p <- sample(1:3, 1)
    ## P1 <- PCG(N)
    Tz <- Toeplitz$new(N)
    case.par <- expand.grid(type = c("fbm", "matern"), b = c(TRUE, FALSE))
    ncase <- nrow(case.par)
    X <- matrix(rnorm(N*p), N, p)
    if(runif(1) < .5) X <- drop(X)
    tol <- 1e-15
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp$type)
      acf <- test_acf_func(N, type)
      Tz$set_acf(acf)
      Tmat <- toeplitz(acf)
      Z <- solve(Tmat, X)
      if(cp$b) {
        expect_equal(Tz$solve(X, method = "pcg", tol = tol), Z)
        expect_equal(solve(Tz, X, method = "pcg", tol = tol), Z)
        ## expect_equal(max(abs((Tmat %*% P1$solve(acf, X, ntol) - X) / Z)), 0, tolerance = 1e-6)
      } else {
        expect_equal(Tmat %*% Tz$solve(method = "pcg", tol = tol), diag(N))
        expect_equal(Tmat %*% solve(Tz, method = "pcg", tol = tol), diag(N))
        ## expect_equal(max(abs((P1$solve(acf, X, ntol)- Z) / Z)), 0, tolerance = 1e-6)
      }
    }
  })
})

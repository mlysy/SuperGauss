source("SuperGauss-testfunctions.R")

context("Toeplitz - Solve.")

nrep <- 10
test_that("GSchur solve method gives correct result.", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 100, sd = 10)))
    d <- sample(1:3, 1)
    Toep <- Toeplitz$new(N)
    case.par <- expand.grid(type = c("fbm", "matern"), b = c(TRUE, FALSE))
    ncase <- nrow(case.par)
    X <- matrix(rnorm(N * d), N, d)
    if(runif(1) > .5) X <- drop(X)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp$type)
      acf <- test_acf_func(N, type)
      Tmat <- toeplitz(acf)
      Z <- solve(Tmat, X)
      Toep$set_acf(acf)
      if(cp$b) {
        expect_equal(solve(Toep, X), Z, tolerance = 1e-5)
        expect_equal(Toep$solve(X), Z, tolerance = 1e-5)
      } else {
        expect_equal(Tmat %*% solve(Toep), diag(N), tolerance = 1e-5)
        expect_equal(Tmat %*% Toep$solve(), diag(N), tolerance = 1e-5)
      }
    }
  })
})

nrep <- 10
test_that("PCG solve method gives correct result.", {
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

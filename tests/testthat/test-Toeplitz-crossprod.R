source("SuperGauss-testfunctions.R")

context("Toeplitz - Cross Products.")

test_that("`t(X) %*% solve(Tz, Y)` is computed correctly.", {
  case_par <- expand.grid(d = 1:3, k = 1:3)
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    cp <- case_par[ii,]
    N <- sample(5:10, 1)
    d <- cp$d
    k <- cp$k
    X <- matrix(rnorm(N * d), N, d)
    Y <- matrix(rnorm(N * k), N, k)
    acf <- exp(-(1:N)/N * runif(1))
    Tz <- toeplitz(acf)
    ldV1 <- 2 * sum(log(diag(chol(Tz))))
    M1 <- crossprod(X, solve(Tz, Y))
    M2 <- DurbinLevinson_crossprod(X, Y, acf, calc_mode = 0)
    expect_equal(M1, M2$IP)
    expect_equal(ldV1, M2$ldV)
  }
})

test_that("`t(X) %*% solve(Tz, X)` is computed correctly.", {
  case_par <- expand.grid(d = 1:3)
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    cp <- case_par[ii,,drop=FALSE]
    N <- sample(5:10, 1)
    d <- cp$d
    ## k <- cp$k
    X <- matrix(rnorm(N * d), N, d)
    ## Y <- matrix(rnorm(N * k), N, k)
    acf <- exp(-(1:N)/N * runif(1))
    Tz <- toeplitz(acf)
    ldV1 <- 2 * sum(log(diag(chol(Tz))))
    M1 <- crossprod(X, solve(Tz, X))
    M2 <- DurbinLevinson_crossprod(X, X, acf, calc_mode = 1)
    expect_equal(M1, M2$IP)
    expect_equal(ldV1, M2$ldV)
  }
})

test_that("`trace(t(X) %*% solve(Tz, Y))` is computed correctly.", {
  case_par <- expand.grid(d = 1:3)
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    cp <- case_par[ii,,drop=FALSE]
    N <- sample(5:10, 1)
    d <- cp$d
    ## k <- cp$k
    X <- matrix(rnorm(N * d), N, d)
    Y <- matrix(rnorm(N * d), N, d)
    acf <- exp(-(1:N)/N * runif(1))
    Tz <- toeplitz(acf)
    ldV1 <- 2 * sum(log(diag(chol(Tz))))
    M1 <- diag(crossprod(X, solve(Tz, Y)))
    M2 <- DurbinLevinson_crossprod(X, Y, acf, calc_mode = 2)
    expect_equal(M1, M2$IP)
    expect_equal(ldV1, M2$ldV)
  }
})

library(SuperGauss)
context("Computation of Positive Definite Toeplitz matrix")

N <- round(abs(rnorm(n = 1, mean = 200, sd = 10)))
d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))

Toep <- new(Toeplitz, N)

test_that("dimension check", {
  expect_identical(Toep$DimCheck(), N)
})

lambda <- rnorm(n = 1, mean = 100, sd = 10)
acf <- exp(-(1:N-1)^2/lambda) * lambda
acf.mat <- toeplitz(acf)
Toep$AcfInput(acf)

X <- matrix(rnorm(N * d), N, d)
zero <- matrix(0, N, d)
Y <- matrix(rnorm(N * d), N, d) * 1e10
Z <- matrix(rnorm(N * d), N, d) * 1e-7

acf0 <- rep(0, N)
acf1 <- rnorm(N)
acf2 <- rnorm(N) * 1e10
acf3 <- rnorm(N) * 1e-7

test_that("Toeplitz-matrix multiplication", {
  expect_equal(Toep$Mult(X), acf.mat %*% X)
  expect_equal(Toep$Mult(zero), acf.mat %*% zero)
  expect_equal(Toep$Mult(Y), acf.mat %*% Y)
  expect_equal(Toep$Mult(Z), acf.mat %*% Z)
})

test_that("Toeplitz-matrix solution", {
  expect_equal(Toep$Solve(X), acf.mat %*% X)
  expect_equal(Toep$Solve(zero), acf.mat %*% zero)
  expect_equal(Toep$Solve(Y), acf.mat %*% Y)
  expect_equal(Toep$Solve(Z), acf.mat %*% Z)
})

test_that("determinant check", {
  expect_equal(Toep$Det(), log(det(acf.mat)))
})

test_that("traceProd check", {
  expect_equal(Toep$TraceProd(acf0), tr(solve(acf.mat, toeplitz(acf0))))
  expect_equal(Toep$TraceProd(acf1), tr(solve(acf.mat, toeplitz(acf1))))
  expect_equal(Toep$TraceProd(acf2), tr(solve(acf.mat, toeplitz(acf2))))
  expect_equal(Toep$TraceProd(acf3), tr(solve(acf.mat, toeplitz(acf3))))
})
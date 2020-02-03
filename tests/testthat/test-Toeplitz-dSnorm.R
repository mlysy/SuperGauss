library(SuperGauss)
library(mvtnorm)
source("test-functions.R")

context("Density")

nrep <- 10
test_that("The GSchur algorithm returns the correct density", {
  replicate(n = nrep, expr = {  
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Mu <- matrix(rnorm(N * d), N, d)
      X <- rSnorm(n = d, acf = acf, fft = FALSE) + Mu
      ld1 <- dSnorm(X = X, mu = Mu, acf = acf, log = TRUE)
      ld2 <- sapply(1:d, function(jj) {
        dmvnorm(x = X[,jj], mean = Mu[,jj], sigma = toeplitz(acf),
                log = TRUE)
      })
      expect_equal(ld1, ld2)
    }
  })
})

test_that("Levinson's algorithm returns the correct density", {
  replicate(n = nrep, expr = {  
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      Mu <- matrix(rnorm(N * d), N, d)
      X <- rSnorm(n = d, acf = acf, fft = FALSE) + Mu
      ld1 <- dSnormDL(X = X, mu = Mu, acf = acf, log = TRUE)
      ld2 <- sapply(1:d, function(jj) {
        dmvnorm(x = X[,jj], mean = Mu[,jj], sigma = toeplitz(acf),
                log = TRUE)
      })
      expect_equal(ld1, ld2)
    }
  })
})

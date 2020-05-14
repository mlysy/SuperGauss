source("SuperGauss-testfunctions.R")

context("Toeplitz - Log-Density (dnormtz wrapper).")

nrep <- 10
test_that("GSchur algorithm returns the correct density", {
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
      X <- rnormtz(n = d, acf = acf, fft = FALSE) + Mu
      ld1 <- dnormtz(X = X, mu = Mu, acf = acf, log = TRUE)
      ld2 <- apply(X-Mu, 2, toep_ldens, gamma = acf)
      ## ld2 <- sapply(1:d, function(jj) {
      ##   mvtnorm::dmvnorm(x = X[,jj], mean = Mu[,jj], sigma = toeplitz(acf),
      ##           log = TRUE)
      ## })
      expect_equal(ld1, ld2)
    }
  })
})

test_that("LTZ algorithm returns the correct density", {
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
      X <- rnormtz(n = d, acf = acf, fft = FALSE) + Mu
      ld1 <- dnormtz(X = X, mu = Mu, acf = acf, log = TRUE, method = "ltz")
      ld2 <- apply(X-Mu, 2, toep_ldens, gamma = acf)
      ## ld2 <- sapply(1:d, function(jj) {
      ##   mvtnorm::dmvnorm(x = X[,jj], mean = Mu[,jj], sigma = toeplitz(acf),
      ##                    log = TRUE)
      ## })
      expect_equal(ld1, ld2)
    }
  })
})

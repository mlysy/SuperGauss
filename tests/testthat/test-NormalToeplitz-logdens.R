source("SuperGauss-testfunctions.R")

context("NormalToeplitz - Log-Density.")

nrep <- 10
test_that("NormalToeplitz$logdens gives correct result.", {
  replicate(n = nrep, expr = {
    N <- round(abs(rnorm(n = 1, mean = 20, sd = 5)))
    case.par <- expand.grid(type = c("fbm", "matern"))
    ncase <- nrow(case.par)
    for(ii in 1:ncase){
      cp <- case.par[ii, ]
      type <- as.character(cp)
      acf <- test_acf_func(N, type)
      X <- rnormtz(n = 1, acf = acf, fft = FALSE)
      Nt <- NormalToeplitz$new(N = N)
      ld1 <- toep_ldens(X, acf)
      ld2 <- Nt$logdens(X, acf)
      expect_equal(ld1, ld2)
    }
  })
})


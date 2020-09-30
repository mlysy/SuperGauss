source("SuperGauss-testfunctions.R")

context("NormalCirculant - Log-Density.")

test_that("NormalCirculant$logdens gives correct result.", {
  skip_on_cran()
  case_par <- expand.grid(N = c(1, sample(2:20, 10)),
                          n_obs = 1:3)
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    list2env(case_par[ii,,drop=FALSE], envir = environment())
    Nu <- floor(N/2)+1
    upsd <- rexp(Nu) * N
    acf <- ifft(unfold_acf(N, upsd))
    uacf <- acf[1:Nu]
    Z <- rnormtz(n = n_obs, acf = acf, fft = FALSE)
    if(N == 1) Z <- t(Z)
    ld1 <- circ_ldens(t(Z), uacf)
    NCt <- NormalCirculant$new(N = N)
    ld2 <- NCt$logdens(Z, uacf)
    expect_equal(ld1, ld2)
  }
})

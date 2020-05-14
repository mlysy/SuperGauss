source("SuperGauss-testfunctions.R")

context("NormalCirculant - Log-Density.")

test_that("NormalCirculant$logdens gives correct result.", {
  case_par <- expand.grid(N = c(1, sample(2:20, 10)))
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    list2env(case_par[ii,,drop=FALSE], envir = environment())
    Nu <- floor(N/2)+1
    upsd <- rexp(Nu) * N
    acf <- ifft(unfold_acf(N, upsd))
    uacf <- acf[1:Nu]
    z <- rnormtz(n = 1, acf = acf, fft = FALSE)
    ld1 <- circ_ldens(z, uacf)
    NCt <- NormalCirculant$new(N = N)
    ld2 <- NCt$logdens(z, uacf)
    expect_equal(ld1, ld2)
  }
})

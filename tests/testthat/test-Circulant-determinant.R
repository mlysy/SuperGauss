source("SuperGauss-testfunctions.R")

context("Circulant - Determiant.")

test_that("Circulant$log_det works as expected.", {
  case_par <- expand.grid(N = c(1, sample(2:20, 10)))
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    list2env(case_par[ii,,drop=FALSE], envir = environment())
    Nu <- floor(N/2)+1
    upsd <- exp(-runif(1) * (1:Nu)/Nu)
    psd <- unfold_acf(N, uacf = upsd)
    acf <- ifft(psd)
    ldet <- sum(log(psd))
    Ct <- Circulant$new(N, uacf = acf[1:Nu])
    expect_equal(ldet, Ct$log_det())
  }
})

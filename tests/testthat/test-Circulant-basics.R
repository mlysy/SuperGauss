source("SuperGauss-testfunctions.R")

context("Circulant - Basics.")

test_that("set_acf/get_acf work as expected.", {
  case_par <- expand.grid(N = c(1, 2, sample(3:20, 10)),
                          sep = c(TRUE, FALSE))
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    list2env(case_par[ii,], envir = environment())
    Nu <- floor(N/2)+1
    uacf <- exp(-runif(1) * (1:Nu)/Nu)
    acf <- unfold_acf(N, uacf = uacf)
    if(sep) {
      # separate set_acf
      Ct <- Circulant$new(N = N)
      Ct$set_acf(uacf)
    } else {
      # combined set_acf
      Ct <- Circulant$new(N = N, uacf = uacf)
    }
    expect_equal(acf, Ct$get_acf())
  }
})

test_that("set_psd/get_psd work as expected.", {
  case_par <- expand.grid(N = c(1, 2, sample(3:20, 10)),
                          input = c("acf", "psd"),
                          sep = c(TRUE, FALSE))
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    list2env(case_par[ii,], envir = environment())
    Nu <- floor(N/2)+1
    upsd <- sort(rexp(Nu), decreasing = TRUE)
    psd <- unfold_acf(N, upsd)
    acf <- ifft(psd)
    uacf <- acf[1:Nu]
    if(input == "acf") {
      if(sep) {
        # separate set_acf
        Ct <- Circulant$new(N = N)
        Ct$set_acf(uacf)
      } else {
        # combined set_acf
        Ct <- Circulant$new(N = N, uacf = uacf)
      }
      expect_equal(psd, Ct$get_psd())
    } else if(input == "psd") {
      if(sep) {
        # separate set_psd
        Ct <- Circulant$new(N = N)
        Ct$set_psd(upsd)
      } else {
        # combined set_psd
        Ct <- Circulant$new(N = N, upsd = upsd)
      }
      expect_equal(acf, Ct$get_acf())
    }
  }
})

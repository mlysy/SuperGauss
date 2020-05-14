source("SuperGauss-testfunctions.R")

context("NormalCirculant - Loglikelihood Gradient (autodiff).")

ntest <- 10
test_that("NormalCirculant$grad_full gives correct result.", {
  case_par <- expand.grid(N = c(1, sample(2:20, ntest-1)))
  ncase <- nrow(case_par)
  for(ii in 1:ncase) {
    list2env(case_par[ii,,drop=FALSE], envir = environment())
    Nu <- floor(N/2)+1
    upsd <- rexp(Nu) * N
    acf <- ifft(unfold_acf(N, upsd))
    uacf <- acf[1:Nu]
    z <- rnormtz(n = 1, acf = acf, fft = FALSE)
    NCt <- NormalCirculant$new(N = N)
    g1_z <- NCt$grad_full(z, uacf, calc_dldz = TRUE, calc_dldu = FALSE)
    g1_u <- NCt$grad_full(z, uacf, calc_dldz = FALSE, calc_dldu = TRUE)
    g1_zu <- NCt$grad_full(z, uacf)
    g2_z <- numDeriv::grad(circ_ldens, x = z, nu = uacf)
    g2_u <- numDeriv::grad(circ_ldens, x = uacf, z = z)
    expect_equal(g2_z, g1_z$dldz)
    expect_equal(g2_z, g1_zu$dldz)
    expect_equal(g2_u, g1_u$dldu)
    expect_equal(g2_u, g1_zu$dldu)
  }
})

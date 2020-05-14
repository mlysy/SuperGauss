source("SuperGauss-testfunctions.R")

context("FFTW - fft/ifft.")

ntest <- 15
test_that("Real-to-Complex fft gives the correct result.", {
  for(ii in 1:ntest) {
    N <- sample(10, 1)
    Nu <- floor(N/2) + 1
    x <- rnorm(N)
    y1 <- fft(x)
    y2 <- SuperGauss:::real_fft(x, inverse = FALSE)
    expect_equal(y1[1:Nu], y2[1:Nu])
  }
})

ntest <- 15
test_that("Complex-to-Real ifft gives the correct result.", {
  for(ii in 1:ntest) {
    N <- sample(10, 1)
    Nu <- floor(N/2) + 1
    y1 <- rnorm(N)
    x <- fft(y1)
    y2 <- SuperGauss:::real_fft(x, inverse = TRUE)
    expect_equal(y1, y2)
  }
})

ntest <- 15
test_that("Real-Even fft gives the correct result.", {
  for(ii in 1:ntest) {
    N <- sample(10, 1)
    Nu <- floor(N/2) + 1
    x <- unfold_acf(N, rnorm(Nu))
    y1 <- Re(fft(x))[1:Nu]
    y2 <- SuperGauss:::even_fft(x, inverse = FALSE)
    expect_equal(y1, y2)
  }
})

ntest <- 15
test_that("Real-Even ifft gives the correct result.", {
  for(ii in 1:ntest) {
    N <- sample(10, 1)
    Nu <- floor(N/2) + 1
    x <- unfold_acf(N, rnorm(Nu))
    y1 <- Re(ifft(x))[1:Nu]
    y2 <- SuperGauss:::even_fft(x, inverse = TRUE)
    expect_equal(y1, y2)
  }
})

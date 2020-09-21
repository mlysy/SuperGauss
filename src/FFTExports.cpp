/// @file FFTExports.cpp
///
/// @brief Wrappers to FFTW functions.

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/RealFFT.h"
#include "SuperGauss/EvenFFT.h"
#include <complex>

/// Compute the FFT of a purely real signal.
///
/// Can also be used to compute the inverse FFT (iFFT) of a complex signal of which the result is a real signal.
///
/// @param[in] x Real vector of length `N` when `inverse = false`.  Otherwise, a complex vector of length `N` having the symmetries required for the iFFT to be purely real.
/// @param[in] inverse Logical flag indicating whether the FFT or iFFT computation is desired.
///
/// @return A complex vector of length `N` when `inverse = false`.  Otherwise a real vector of length `N`.
///
// [[Rcpp::export]]
SEXP real_fft(SEXP x, bool inverse = false) {
  if(!inverse) {
    NumericVector x_(x);
    int N = x_.length();
    std::complex<double>* y_ = new std::complex<double>[N];
    ComplexVector y(N);
    RealFFT rfft(N);
    rfft.fft(y_, REAL(x_));
    for(int ii=0; ii<N; ii++) {
      y[ii].r = y_[ii].real();
      y[ii].i = y_[ii].imag();
    }
    return y;
  } else {
    ComplexVector xcv_(x);
    int N = xcv_.length();
    std::complex<double>* x_ = new std::complex<double>[N];
    NumericVector y(N);
    RealFFT rfft(N);
    for(int ii=0; ii<N; ii++) {
      x_[ii] = std::complex<double>(xcv_[ii].r, xcv_[ii].i);
    }
    rfft.ifft(REAL(y), x_);
    return y;
  }
}

/// Compute the FFT of an even signal.
///
/// Zero-centered even vectors are real vectors `x = (x_0, ..., x_{N-1})` of the form
///
/// ```
/// x = (x_0, ..., x_{Nu-1}, x_{Nu-2}, ..., x_1),   N = 2*Nu-2   even,
///   = (x_0, ..., x_{Nu-1}, x_{Nu-1}, ..., x_1),   N = 2*Nu-1   odd,
/// ```
///
/// where `Nu = floor(n/2) + 1`.  The FFT and iFFT of `x` are also zero-centered even vectors and thus benefit from additional optimization.
///
/// @param[in] x Even vector of length `N`.  Only the first `Nu` elements will be used.
/// @param[in] inverse Logical; whether to compute the FFT or its inverse.
///
/// @return The FFT or iFFT of `x`.  Only the first `Nu` elements are returned.
///
// [[Rcpp::export]]
NumericVector even_fft(NumericVector x, bool inverse = false) {
  int N = x.length();
  EvenFFT efft(N);
  int Nu = efft.usize(); // storage size
  NumericVector y(Nu);
  if(!inverse) {
    efft.fft(REAL(y), REAL(x));
  } else {
    efft.ifft(REAL(y), REAL(x));
  }
  return y;
}


/// @file FFTExports.cpp

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/RealFFT.h"
#include "SuperGauss/EvenFFT.h"
#include <complex>

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


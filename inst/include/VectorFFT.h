/// @file VectorFFT.h

///////////////////////////////////////////////
// Fast Fourier Transformation and Inverse
///////////////////////////////////////////////

#ifndef VectorFFT_h
#define VectorFFT_h 1

// usual header
#include <complex>
#include <algorithm>
#include <fftw3.h>
// #include <Rcpp.h>
// #include <iostream>

/// @brief FFT for real to complex and corresponding iFFT for complex to real.
///
/// Allocates memory for the corresponding `fftw` operations within the object, copies memory in and out each time the FFT and iFFT members are called.
class VectorFFT {
private:
  typedef std::complex<double> dcomplex; ///< Typedef for complex double
  fftw_plan planfwd_;  ///< Plan for forward FFT.
  fftw_plan planback_;  ///< Plan for backward FFT.
  fftw_complex* y_; ///< Where to compute FFT.
  double* x_; ///< Where to compute iFFT.
  int n_; ///< Size of input vector.
  int npad_; ///< Padded size of input.
public:
  /// Constructor.
  VectorFFT(int n);
  /// Destructor.
  ~VectorFFT();
  /// Perform the FFT on the input data.
  void fft(std::complex<double>* y, const double* x);
  /// Perform the inverse FFT on the input data.
  void ifft(double* x, const std::complex<double>* y);
  /// Get size of FFT.
  int size() {return n_;}
};

/// @param[in] n Size of FFT/iFFT to be computed.
inline VectorFFT::VectorFFT(int n) {
  n_ = n;
  npad_ = ceil((double)(n + 1) / 2);
  // x_ = new double[n];
  x_ = fftw_alloc_real(n);
  std::fill(x_, x_ + n, 0);
  // y_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  y_ = fftw_alloc_complex(n);
  planfwd_ = fftw_plan_dft_r2c_1d(n, x_, y_, FFTW_ESTIMATE);
  planback_ = fftw_plan_dft_c2r_1d(n, y_, x_, FFTW_ESTIMATE);
  return;
}

inline VectorFFT::~VectorFFT() {
  // delete[] x_;
  fftw_free(x_);
  fftw_free(y_);
  fftw_destroy_plan(planfwd_);
  fftw_destroy_plan(planback_);
}

/// Calculates `y = FFT(x)`.
///
/// @param[out] y Complex vector output.
/// @param[in] x  Real vector input.
inline void VectorFFT::fft(std::complex<double>* y,
			   const double* x) {
  std::copy(x, x + n_, x_);
  fftw_execute(planfwd_);
  for (int ii = 0; ii < npad_; ++ii) {
    y[ii] = dcomplex(y_[ii][0], y_[ii][1]);
  }
  return;
}

/// Calculates `x = iFFT(y)`.
///
/// @param[out] x Real vector output.
/// @param[in] y Complex vector input.
inline void VectorFFT::ifft(double* x, const std::complex<double>* y) {
  for (int ii = 0; ii < npad_; ++ii) {
    y_[ii][0] = real(y[ii]);
    y_[ii][1] = imag(y[ii]);
  }

  fftw_execute(planback_);
  for (int ii = 0; ii < n_; ++ii) {
    x[ii] = x_[ii] / n_;
  }
  return;
}

#endif

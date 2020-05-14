/// @file RealFFT.h

#ifndef RealFFT_h
#define RealFFT_h 1

// usual header
#include <complex>
#include <algorithm>
#include <fftw3.h>
// #include <Rcpp.h>
// #include <iostream>

/// @brief FFT for real to complex and corresponding iFFT for complex to real.
///
/// Allocates memory for the corresponding `fftw` operations within the object, copies memory in and out each time the FFT and iFFT members are called.
///
/// The FFT of a real vector `x = (x_0, ..., x_{N-1})` produces a complex vector `y = (y_0, ..., y_{N-1})` such that `y_k = conj(y_{k mod N})`.  Therefore, `y` is completely determined from its first `Nu = floor(N/2)+1` elements, such that only these `Nu` elements are returned as FFT output and queried as iFFT input.
class RealFFT {
private:
  typedef std::complex<double> dcomplex; ///< Typedef for complex double
  fftw_plan planfwd_;  ///< Plan for forward FFT.
  fftw_plan planback_;  ///< Plan for backward FFT.
  fftw_complex* y_; ///< Where to compute FFT.
  double* x_; ///< Where to compute iFFT.
  int N_; ///< Size of input vector.
  int Nu_; ///< Padded size of input.
public:
  /// Constructor.
  RealFFT(int N);
  /// Destructor.
  ~RealFFT();
  /// Perform the FFT on the input data.
  void fft(std::complex<double>* y, const double* x);
  /// Perform the inverse FFT on the input data.
  void ifft(double* x, const std::complex<double>* y);
  /// Get size of FFT.
  int size() {return N_;}
};

/// @param[in] N Size of FFT/iFFT to be computed.
inline RealFFT::RealFFT(int N) {
  N_ = N;
  Nu_ = ceil((double)(N_ + 1) / 2);
  // x_ = new double[n];
  x_ = fftw_alloc_real(N_);
  std::fill(x_, x_ + N_, 0);
  // y_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  y_ = fftw_alloc_complex(N_);
  planfwd_ = fftw_plan_dft_r2c_1d(N_, x_, y_, FFTW_ESTIMATE);
  planback_ = fftw_plan_dft_c2r_1d(N_, y_, x_, FFTW_ESTIMATE);
}

inline RealFFT::~RealFFT() {
  // delete[] x_;
  fftw_free(x_);
  fftw_free(y_);
  fftw_destroy_plan(planfwd_);
  fftw_destroy_plan(planback_);
}

/// Calculates `y = FFT(x)`.
///
/// @param[out] y Complex vector output of length `Nu = floor(N/2)+1`.
/// @param[in] x  Real vector input of length `N`.
inline void RealFFT::fft(std::complex<double>* y,
			   const double* x) {
  std::copy(x, x + N_, x_);
  fftw_execute(planfwd_);
  for (int ii = 0; ii < Nu_; ++ii) {
    y[ii] = dcomplex(y_[ii][0], y_[ii][1]);
  }
  return;
}

/// Calculates `x = iFFT(y)`.
///
/// @param[out] x Real vector output of length `N`.
/// @param[in] y Complex vector input of length `Nu = floor(N/2)+1`.
inline void RealFFT::ifft(double* x, const std::complex<double>* y) {
  for (int ii = 0; ii < Nu_; ++ii) {
    y_[ii][0] = y[ii].real();
    y_[ii][1] = y[ii].imag();
  }

  fftw_execute(planback_);
  for (int ii = 0; ii < N_; ++ii) {
    x[ii] = x_[ii] / N_;
  }
  return;
}

#endif

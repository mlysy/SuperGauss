/// @file EvenFFT.h

/// FFT and iFFT for zero-centered even vectors.

#ifndef EvenFFT_h
#define EvenFFT_h 1

// usual header
#include <complex>
#include <algorithm>
#include <fftw3.h>

/// @brief FFT and iFFT for zero-centered even vectors.
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
/// This class allocates memory for the corresponding `fftw` operations within the object, and copies memory in and out each time the FFT and iFFT members are called.
class EvenFFT {
 private:
  fftw_plan planeven_;  ///< Plan for even FFT/iFFT.
  fftw_plan planodd_;  ///< Plan for odd FFT/iFFT.
  double* x_; ///< Input for FFT/iFFT.
  double* yeven_; ///< Output for even FFT/iFFT.
  fftw_complex* yodd_; ///< Output for odd FFT/iFFT.
  int N_; ///< Size of input/output vectors.
  int Nu_; ///< Number of unique input/output elements.
  bool Neven_; ///< Whether the input/output size is an even integer.
public:
  /// Constructor.
  EvenFFT(int N);
  /// Destructor.
  ~EvenFFT();
  /// Perform the FFT on the input data.
  void fft(double* y, const double* x);
  /// Perform the inverse FFT on the input data.
  void ifft(double* x, const double* y);
  /// Get size of input/output.
  int size() {
    return N_;
  }
  /// Get number of unique elements of input/outpu.
  int usize() {
    return Nu_;
  }
};

/// @param[in] N Size of FFT/iFFT to be computed.
inline EvenFFT::EvenFFT(int N) {
  N_ = N;
  Nu_ = N_/2 + 1;
  Neven_ = (N%2 == 0);
  if(Neven_) {
    x_ = fftw_alloc_real(Nu_);
    yeven_ = fftw_alloc_real(Nu_);
    planeven_ = fftw_plan_r2r_1d(Nu_, x_, yeven_, FFTW_REDFT00, FFTW_ESTIMATE);
  } else {
    x_ = fftw_alloc_real(N_);
    yodd_ = fftw_alloc_complex(N_);
    planodd_ = fftw_plan_dft_r2c_1d(N_, x_, yodd_, FFTW_ESTIMATE);
  }
}

inline EvenFFT::~EvenFFT() {
  fftw_free(x_);
  if(Neven_) {
    fftw_free(yeven_);
    fftw_destroy_plan(planeven_);
  } else {
    fftw_free(yodd_);
    fftw_destroy_plan(planodd_);
  }
}

/// Calculates `y = FFT(x)`.
///
/// @param[out] y Real vector output.
/// @param[in] x Real vector input.
inline void EvenFFT::fft(double* y, const double* x) {
  std::copy(x, x + Nu_, x_);
  if(Neven_) {
    fftw_execute(planeven_);
    std::copy(yeven_, yeven_ + Nu_, y);
  } else {
    std::reverse_copy(x + 1, x + Nu_, x_ + Nu_);
    fftw_execute(planodd_);
    for(int ii=0; ii<Nu_; ii++) {
      y[ii] = yodd_[ii][0];
    }
  }
  return;
}

/// Calculates `x = iFFT(y)`.
///
/// @param[out] x Real vector output.
/// @param[in] y Real vector input.
inline void EvenFFT::ifft(double* x, const double* y) {
  fft(x, y);
  for(int ii=0; ii<Nu_; ii++) {
    x[ii] /= double(N_);
  }
  // std::copy(y, y + Nu_, y_);
  // fftw_execute(planback_);
  // for (int ii = 0; ii < npad_; ++ii) {
  //   x[ii] = x_[ii] / n_;
  // }
  return;
}


#endif

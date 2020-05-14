/// @file Circulant.h

#ifndef Circulant_h
#define Circulant_h 1

#include "RealFFT.h"
#include "EvenFFT.h"
#include <complex>
#include <algorithm>

/// @brief Circulant matrix class.
///
/// Provides linear algebra methods for Circulant symmetric positive-definite (i.e., variance) matrices.  An `N x N` circulant variance matrix `Ct` is defined by the first `Nu = floor(N/2)+1` unique elements of its first row/column (i.e., autocorrelation vector).  That is, the first row of `Ct` is given by
///
/// ```
/// Ct[0,:] = (uacf_0, ..., uacf_{Nu-1}, uacf_{Nu-2}, ..., uacf_1),  N = 2*Nu-2,
///         = (uacf_0, ..., uacf_{Nu-1}, uacf_{Nu-1}, ..., uacf_1),  N = 2*Nu-1.
/// ```
///
/// @note No check is made as to whether `Ct` is positive-definite. Non-negative matrix power product works for non positive-definite `Ct`, but the behavior of log-determinant in that case is undefined.
class Circulant {
 private:
  typedef std::complex<double> dcomplex;
  // internal storage.  don't overwrite locally.
  int N_; ///< Size of Circulant matrix.
  int Nu_; ///< Number of unique elements.
  bool Neven_; ///< Whether N is even.
  double* acf_; ///< Autocorrelation vector of the Circulant matrix.
  double* psd_; ///< FFT of the autocorrelation vector.
  double logdet_; ///< Storage for log-determinant of Circulant matrix.
  EvenFFT* efft_; ///< Object for acf FFT computations.
  RealFFT* rfft_; ///< Object for real-valued FFT computations.
  bool has_acf_;    ///< Wheter input argument `acf_` has been modified.
  bool has_psd_;   ///< Whether the PSD has been computed.
  bool has_logdet_;   ///< Whether log-determinant has been computed.
  dcomplex* x_fft_; // temporary FFT storage
  /// Unfold a vector.
  void unfold(double* x);
  /// Copy and unfold a vector.
  void unfold_copy(double* x, const double* xu);
  /// Precomputations for matrix-vector multiplication.
  void prod_setup();
public:
  /// Constructor.
  Circulant(int N);
  /// Destructor.
  ~Circulant();
  /// Set the acf of the Circulant matrix.
  void set_acf(const double* uacf);
  /// Get the acf of the Circulant matrix.
  void get_acf(double* acf);
  /// Set the psd of the Circulant matrix.
  void set_psd(const double* upsd);
  /// Void get the psd of the Circulant matrix.
  void get_psd(double* psd);
  /// Size of the Circulant matrix.
  int size();
  /// Number of unique elements in the Circulant matrix.
  int usize();
  /// Check whether the acf has been set.
  bool has_acf();
  /// Circulant matrix-vector multiplication.
  void prod(double* y, const double* x);
  /// Circulant matrix-vector power product.
  void prod_pow(double* y, const double* x, double alpha);
  /// Solve Circulant matrix-vector system of equations.
  void solve(double* y, const double* x);
  /// Log-determinant of the Circulant matrix.
  double log_det();
};

/// Copies the first half of x to the second half in reverse order.
///
/// @param[in/out] x Vector of which `Nu` elements are input and `N` elements are output
inline void Circulant::unfold(double* x) {
  std::reverse_copy(x+1, x + Nu_ - Neven_, x + Nu_);
  return;
}


inline void Circulant::unfold_copy(double* x, const double* xu) {
  std::copy(xu, xu + Nu_, x); // first half
  unfold(x); // second half
  // std::reverse_copy(xu+1, xu + Nu_ - Neven_, x + Nu_); // second half
  return;
}

/// @param[in] N Size of Circulant matrix.
inline Circulant::Circulant(int N) {
  N_ = N;
  Nu_ = N_/2 + 1;
  Neven_ = (N_%2 == 0);
  acf_ = new double[N_];
  psd_ = new double[N_];
  x_fft_ = new dcomplex[N_];
  rfft_ = new RealFFT(N_);
  efft_ = new EvenFFT(N_);
  has_acf_ = false;
  has_psd_ = false;
  has_logdet_ = false;
}

inline Circulant::~Circulant() {
  delete [] acf_;
  delete[] psd_;
  delete[] x_fft_;
  delete efft_;
  delete rfft_;
}

/// @return The size `N` of the Circulant matrix.
inline int Circulant::size() {
  return N_;
}

/// @return The number of unique elements `Nu = floor(N/2)+1`.
inline int Circulant::usize() {
  return Nu_;
}

/// @param[in] uacf Unique elements of first row/column of Circulant matrix.
inline void Circulant::set_acf(const double* uacf) {
  unfold_copy(acf_, uacf);
  // std::copy(uacf, uacf + Nu_, acf_); // first half
  // std::reverse_copy(uacf+1, uacf + Nu_ - Neven_, acf_ + Nu_); // second half
  has_acf_ = true;
  has_psd_ = false;
  has_logdet_ = false;
  return;
}

/// @param[out] acf Complete autocorrelation vector of length `N`.
inline void Circulant::get_acf(double* acf) {
  std::copy(acf_, acf_ + N_, acf);
  return;
}

/// The power spectral density (PSD) of the autocorrelation vector is defined as `psd = FFT(acf)`.  For Circulant matrices, both `psd` and `acf` are nonnegative real vectors with even symmetry. 
///
/// @param[in] upsd The `Nu` unique elements of the PSD.
inline void Circulant::set_psd(const double* upsd) {
  // // std::copy(upsd, upsd + Nu_, psd_); // first half
  // // std::reverse_copy(upsd+1, upsd + Nu_ - Neven_, psd_ + Nu_); // second half
  // for(int ii=0; ii<Nu_; ii++) {
  //   psd_[ii].real(upsd[ii]);
  //   psd_[ii].imag(0.0);
  // }
  // for(int ii=1; ii<Nu_ - Neven_; ii++) {
  //   psd_[N_-ii].real(upsd[ii]);
  //   psd_[N_-ii].imag(0.0);
  // }
  unfold_copy(psd_, upsd);
  // convert psd to acf
  efft_->ifft(acf_, psd_);
  unfold(acf_);
  has_acf_ = true;
  has_psd_ = true;
  has_logdet_ = false;
  return;
}

/// @param[out] Complete PSD vector of length `N`.
inline void Circulant::get_psd(double* psd) {
  if(!has_psd_) prod_setup();
  // for(int ii=0; ii<N_; ii++) {
  //   psd[ii] = psd_[ii].real();
  // }
  std::copy(psd_, psd_ + N_, psd);
  return;
}

/// @return `true` if `set_acf()` has been called, `false` otherwise.
inline bool Circulant::has_acf() {
  return has_acf_;
}

/// Precomputes the FFT of `acf_` into `psd_` and sets `has_psd_ = true`.
inline void Circulant::prod_setup() {
  // printf("prod_setup\n");
  // calculate FFT
  efft_->fft(psd_, acf_);
  unfold(psd_);
  has_psd_ = true;
  // for(int ii=0; ii<N_; ii++) {
    // printf("acf_[%i] = %f\n", ii, acf_[ii]);
    // printf("psd_[%i] = (%f, %f)\n", ii, psd_[ii].real(), psd_[ii].imag());
    // acf_ is even so psd_ is purely real
    // psd_[ii][1] = 0.0; 
  //   psd_[ii].imag(0.0);
  // }
  return;
}

/// Computes `y = Ct * x`, where `Ct = Circulant(uacf)`.  This is done in `O(N log N)` operations using the FFT.
///
/// @param[out] y Output vector of size `N` for the matrix-vector multiplication `Ct * x`.
/// @param[in] x Input vector of size `N`.
inline void Circulant::prod(double* y, const double* x) {
  if(!has_psd_) prod_setup();
  rfft_->fft(x_fft_, x);
  for(int ii=0; ii<N_; ii++) {
    x_fft_[ii] *= psd_[ii];
    // printf("x_fft_[%i] = (%f, %f)\n", ii, x_fft_[ii].real(), x_fft_[ii].imag());
  }
  rfft_->ifft(y, x_fft_);
  // for(int ii=0; ii<N_; ii++) printf("y[%i] = %f\n", ii, y[ii]);
  return;
}

/// Computes `y = Ct^{-1} x`, where `Ct = Circulant(uacf)`.  Can be solved in `O(N log N)` steps using the FFT.
///
/// @param[out] y Vector of length `N` containing the output `y = Ct^{-1} x`.
/// @param[in] x Input vector of length `N`.
inline void Circulant::solve(double* y, const double* x) {
  if(!has_psd_) prod_setup();
  rfft_->fft(x_fft_, x);
  for(int ii=0; ii<N_; ii++) {
    x_fft_[ii] /= psd_[ii];
  }
  rfft_->ifft(y, x_fft_);
  return;
}

/// Computes `y = Ct^{alpha} * x`, where `Ct = Circulant(uacf)`.  This is done in `O(N log N)` operations using the FFT.
///
/// @param[out] y Vector of size `N` containing the output `y = Ct^{alpha} x`.
/// @param[in] x Input vector of size `N`.
/// @param[in] alpha Real valued matrix exponent.
inline void Circulant::prod_pow(double* y, const double* x, double alpha) {
  if(!has_psd_) prod_setup();
  rfft_->fft(x_fft_, x);
  for(int ii=0; ii<N_; ii++) {
    x_fft_[ii] *= std::pow(psd_[ii], alpha);
  }
  rfft_->ifft(y, x_fft_);
  return;
}

/// @return Scalar containing the value of `log(det(Circulant(uacf)))`.
inline double Circulant::log_det() {
  if(!has_psd_) prod_setup();
  if(!has_logdet_) {
    logdet_ = 0.0;
    for(int ii=0; ii<N_; ii++) {
      logdet_ += std::log(psd_[ii]);
    }
    has_logdet_ = true;
  }
  return logdet_;
}

#endif

/// @file Toeplitz.h

///////////////////////////////////////////////
// Toeplitz matrix class
///////////////////////////////////////////////

#ifndef Toeplitz_h
#define Toeplitz_h 1

#include "GSchur.h"


/// @brief Toeplitz matrix class.
///
/// Provides linear algebra methods for Toeplitz symmetric positive-definite (i.e., variance) matrices including matrix multiplication, inversion, determinant, and first and second derivatives of the log-determinant.
class Toeplitz {
private:
  typedef std::complex<double> dcomplex;
  // internal storage.  don't overwrite locally.
  int N_;        ///< Size of Toeplitz matrix.
  int N2_;       ///< Size of FFT product.
  double* acf_;    ///< First column of the Toeplitz matrix.
  double* tzcirc_; ///< Storage for Toeplitz circulant embedding.
  dcomplex* tzcirc_fft_; ///< FFT of Toeplitz circulant embedding.
  double* delta_; ///< Storage for first column of the inverse-Toeplitz matrix.
  double logdet_; ///< Storage for log-determinant of Toeplitz matrix.
  double traceinv_;///< Trace of inverse-Toeplitz.
  dcomplex* conv_fft_; ///< FFT for convolutions.
  // double* vec_zero_; ///< Storage for zero-padding
  GSchurN* gs_;     ///< Object to compute GSchur algorithm.
  RealFFT* rfft_;   ///< Object for fft computations.
  bool has_acf_;    ///< Wheter input argument `acf_` has been modified.
  bool has_prod_;   ///< Whether multiplication-related FFT has been done.
  bool has_solve_;  ///< Whether inversion-related FFT has been done.
  bool has_trace_;  ///< Whether trace-of-inverse calculation has been done.
  dcomplex* L1_fft_; ///< FFT for first lower-triangular Toeplitz matrix in `solve`.
  dcomplex* tL1_fft_; ///< FFT for transpose of first lower-triangular Toeplitz matrix in `solve`.
  dcomplex* L2_fft_; ///< FFT for second lower-triangular Toeplitz matrix in `solve`.
  dcomplex* tL2_fft_; ///< FFT for transpose of second lower-triangular Toeplitz matrix in `solve`.
  // temporary storage.  ok to overwrite locally.
  double *vec1_, *vec2_, *vec3_, *vec4_, *vec5_, *vec6_;
  dcomplex *vec1_fft_, *vec2_fft_, *vec3_fft_; //, *vec4_fft_, *vec5_fft_;
  /// Precomputations for matrix-vector multiplication.
  void prod_setup();    
  /// Precomputations for solving linear systems.
  void solve_setup();
  /// Trace of product of lower and upper triangular Toeplitz matrices.
  double trace_LU(const double* L, const double* U);
  /// Convolution of real vectors from their Fourier inputs.
  void conv_fft(double* z, const dcomplex* x_fft, const dcomplex* y_fft);
  /// Zero-padded FFT of real vector.
  void zero_fft(dcomplex* y_fft, double* x);
public:
  /// Constructor.
  Toeplitz(int N, int bmod);  
  /// Destructor.
  ~Toeplitz();
  /// Set the acf of the Toeplitz matrix.
  void set_acf(const double* acf);
  /// Get the acf of the Toepliz matrix.
  void get_acf(double* acf);  
  /// Size of the Toeplitz matrix.
  int size(); 
  /// Check whether the acf has been set.
  bool has_acf();  
  /// Toeplitz matrix-vector multiplication.
  void prod(double* y, const double* x);
  /// External symmetric Toeplitz matrix-vector multiplication.
  void prod(double* y, const double* x, const double* acf1);
  /// External non-symmetric Toeplitz matrix-vector multiplication.
  void prod(double* y, const double* x, const double* acf1, const double* acf2);  
  /// Solve Toeplitz matrix-vector system of equations.
  void solve(double* y, const double* x);
  /// Log-determinant of the Toeplitz matrix.
  double log_det();
  /// Trace of inverse-Toeplitz matrix.
  double trace_inv();
  /// Gradient-specialized trace-product.
  double trace_grad(const double* acf1);
  /// Hessian-specialized trace-product.
  double trace_hess(const double* acf1, const double* acf2);
};

/// @param[in] N Size of Toeplitz matrix.
/// @param[in] bmod Size of binary modulus for GSchur calculation.
inline Toeplitz::Toeplitz(int N, int bmod = 64) {
  N_ = N;
  N2_ = 2 * (N_ / 2 + 1);
  // N3_ = N_ + 1;
  acf_ = new double[N_];
  has_acf_ = false;
  has_prod_ = false;
  has_solve_ = false;
  has_trace_ = false;
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1
    gs_ = new GSchurN(N_, bmod);
    rfft_ = new RealFFT(2 * N_);
    tzcirc_ = new double[2 * N_];
    tzcirc_fft_ = new dcomplex[2 * N_];
    conv_fft_ = new dcomplex[2 * N_];
    // vec_zero_ = new double[2 * N_];
    delta_ = new double[N_];
    L1_fft_ = new dcomplex[2 * N_];
    tL1_fft_ = new dcomplex[2 * N_];
    L2_fft_ = new dcomplex[2 * N_];
    tL2_fft_ = new dcomplex[2 * N_];
    vec1_ = new double[2 * N_];
    vec1_fft_ = new dcomplex[2 * N_];
    vec2_ = new double[2 * N_];
    vec2_fft_ = new dcomplex[2 * N_];
    vec3_ = new double[2 * N_];
    vec3_fft_ = new dcomplex[2 * N_];
    vec4_ = new double[2 * N_];
    // vec4_fft_ = new dcomplex[2 * N_];
    vec5_ = new double[2 * N_];
    // vec5_fft_ = new dcomplex[2 * N_];
    vec6_ = new double[2 * N_];
  }
}

inline Toeplitz::~Toeplitz() {
  delete[] acf_;
  // GSchur algorithm only supports N > 1 case.
  if (N_ > 1) {
    delete gs_;
    delete rfft_;
    delete[] tzcirc_;
    delete[] tzcirc_fft_;
    delete[] conv_fft_;
    delete[] delta_;
    // delete[] vec_zero_;
    delete[] L1_fft_;
    delete[] tL1_fft_;
    delete[] L2_fft_;
    delete[] tL2_fft_;
    delete[] vec1_;
    delete[] vec1_fft_;
    delete[] vec2_;
    delete[] vec2_fft_;
    delete[] vec3_;
    delete[] vec3_fft_;
    delete[] vec4_;
    // delete[] vec4_fft_;
    delete[] vec5_;
    // delete[] vec5_fft_;
    delete[] vec6_;    
  }
}

/// @param[in] acf First row/column of Toeplitz matrix.
///
/// @note Calls to `Toeplitz::prod()`, `Toeplitz::solve()`, and `Toeplitz::trace_inv()` store intermediate calculations which make these calls much faster for the same `acf` with different values of the other inputs.  Calling `set_acf()` indicates to `Toeplitz` that these intermediate calculations need to be recomputed.
inline void Toeplitz::set_acf(const double* acf) {
  std::copy(acf, acf + N_, acf_);
  has_acf_ = true;
  has_prod_ = false;
  has_solve_ = false;
  has_trace_ = false;
  return;
}

/// @param[out] acf First row/column of Toeplitz matrix.
inline void Toeplitz::get_acf(double* acf) {
  std::copy(acf_, acf_ + N_, acf);
  return;
}

inline bool Toeplitz::has_acf() { 
  return has_acf_; 
}

inline int Toeplitz::size() { 
  return N_; 
}

/// Precomputes the FFT of the circulant embedding of `acf_` into `tzcirc_fft_` and sets `has_prod_ = true`.
inline void Toeplitz::prod_setup() {
  has_prod_ = true;
  if (N_ > 1) {
    std::copy(acf_, acf_ + N_, tzcirc_);
    tzcirc_[N_] = 0;
    std::reverse_copy(acf_ + 1, acf_ + N_, tzcirc_ + N_ + 1);
    // std::copy(acf_ + 1, acf_ + N_, tzcirc_ + N_ + 1);
    // std::reverse(tzcirc_ + N_ + 1, tzcirc_ + 2 * N_);
    rfft_->fft(tzcirc_fft_, tzcirc_);
  }
  return;
}

/// Calculates the convolution between length-`N_` vectors `x` and `y` based on their (zero-padded) FFT transformations `x_fft` and `y_fft`.  This is an elementwise multiplication of `x_fft` and `y_fft` followed by an iFFT back to the real domain.
///
/// @param[out] z Convolution output.
/// @param[in] x_fft FFT of first input vector.
/// @param[in] y_fft FFT of second input vector.
inline void Toeplitz::conv_fft(double* z, const dcomplex* x_fft,
			       const dcomplex* y_fft) {
  for(int ii=0; ii<N2_; ii++) {
    conv_fft_[ii] = x_fft[ii] * y_fft[ii];
  }
  // complex_mult(conv_fft_, x_fft, y_fft, N2_);
  rfft_->ifft(z, conv_fft_);
  return;
}

/// Calculates the FFT of the length-`2N_` vector `x_zero = [x, 0, ..., 0]`.
///
/// @param[out] y_fft FFT of zero-padded vector `x_zero`.
/// @param[in] x Input vector.
///
/// @warning The input vector `x` is modified such that `x[i] = 0` for `i=N,...,2N-1`.
inline void Toeplitz::zero_fft(dcomplex* y_fft, double* x) {
  std::fill(x + N_, x + 2 * N_, 0.0);
  rfft_->fft(y_fft, x);
}

/// Calculates `trace(L * U)`, where `L` and `U` are lower/upper triangular Toeplitz matrices. It is an `O(N)` algorithm.
///
/// @param[in] L First column of `L`.
/// @param[in] U First row of `U`.
/// @return The trace-product `trace(L * U)`.
inline double Toeplitz::trace_LU(const double* L, const double* U) {
  double trace = 0;
  for (int ii = 0; ii < N_; ++ii) {
    trace += (N_ - ii) * L[ii] * U[ii];
  }
  return trace;
}


/// Computes `y = Tz * x`, where `Tz = Toeplitz(acf)`.  This is done in `O(N log N)` operations by embedding `Tz` into a circulant matrix, zero-padding `x`, and performing the circulant matrix-vector multiplication using the FFT.
///
/// @param[out] y Output vector of size `N` for the matrix-vector multiplication `Tz * x`.
/// @param[in] x Input vector of size `N`.
inline void Toeplitz::prod(double* y, const double* x) {
  // Pointers to temporary storage: x_, x_fft_, y_.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  // dcomplex* y_fft_ = vec2_fft_;
  // tzcirc_ = [acf_, 0, rev(acf_[-1])]
  if(!has_prod_) prod_setup();
  std::copy(x, x + N_, x_);
  zero_fft(x_fft_, x_);
  // y_ = ifft(fft(tzcirc_) * fft(x_))[1:N_]
  conv_fft(y_, tzcirc_fft_, x_fft_);
  std::copy(y_, y_ + N_, y);
  return;
}


/// Precomputes the FFT of the four triangular-Toeplitz matrices involved in the Gohberg-Semencul formula, and sets `has_solve_ = true`.
inline void Toeplitz::solve_setup() {
  // Pointers to temporary storage: z_.
  double* z_ = vec1_;
  has_solve_ = true;
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    gs_->compute(delta_, logdet_, acf_);
    // tL1_fft_ stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_1'
    z_[0] = delta_[0];
    std::fill(z_ + 1, z_ + N_ + 1, 0.0);
    std::reverse_copy(delta_ + 1, delta_ + N_, z_ + N_ + 1);
    rfft_->fft(tL1_fft_, z_);
    // L1_fft_ stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_1
    std::copy(delta_, delta_ + N_, z_);
    zero_fft(L1_fft_, z_);
    // tL2_fft_ stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_2'
    std::fill(z_, z_ + N_ + 1, 0.0);
    std::copy(delta_ + 1, delta_ + N_, z_ + N_ + 1);
    rfft_->fft(tL2_fft_, z_);	
    // L2_fft_ stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_2
    std::fill(z_, z_ + 2 * N_, 0.0);
    std::reverse_copy(delta_ + 1, delta_ + N_, z_ + 1);
    rfft_->fft(L2_fft_, z_);
  }
  return;
}

/// Linear systems with symmetric positive-definite Toeplitz matrices `Tz = Toeplitz(acf)` can be solved in `O(N log^2 N)` steps using the approach of Ammar & Gragg (1988):
///
/// 1.  Use the Generalized Schur (GSchur) algorithm to obtain the first column `delta` of `Tz^{-1}` in `O(N log^2 N)` steps.
/// 2.  The first column `delta` of `Tz^{-1}` can be used to obtain the Gohberg-Semencul decomposition
///    ```
///    Tz^{-1} = 1/rho (L1 L1' - L2 L2'),
///    ```
///    where the scalar `rho` and the lower-triangular Toeplitz matrices `L1` and `L2` are all functions of `delta`.  Using this decomposition, the calculation of `y = Tz^{-1} x` involves only Toeplitz matrix-vector multiplications, each of which is `O(N log N)` with the FFT circulant embedding.
///
/// @param[out] y Vector of length `N` containing the output `y = Tz^{-1} x`.
/// @param[in] x Input vector of length `N`.
inline void Toeplitz::solve(double* y, const double* x) {
  // Pointers to temporary storage: x_, x_fft_, y_, y_fft_, z_.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  dcomplex* y_fft_ = vec2_fft_;
  double* z_ = vec3_;
  // dcomplex* z_fft_ = vec3_fft_;
  if (!has_solve_) solve_setup();
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(x, x + N_, x_);
    zero_fft(x_fft_, x_);
    // y_ = L1' * x_
    conv_fft(y_, tL1_fft_, x_fft_);
    // z_ = L1 * y_ = L1 * L1' * x_
    zero_fft(y_fft_, y_);
    conv_fft(z_, L1_fft_, y_fft_);
    // y_ = L2' * x_
    conv_fft(y_, tL2_fft_, x_fft_);
    // x_ = L2 * y_ = L2 * L2' * x_ (temporarily stored in x_)
    zero_fft(y_fft_, y_);
    conv_fft(x_, L2_fft_, y_fft_);
    // y = (z_ - x_) / delta[1] = 1/delta[1] * (L1 * L1' * x_ - L2 * L2' * x_)
    for (int ii = 0; ii < N_; ++ii) {
      y[ii] = (z_[ii] - x_[ii]) / delta_[0];
    }
  } else {
    // N = 1 case.
    y[0] = x[0] / acf_[0];
  }
  return;
}

/// The log-determinant of `Tz = Toeplitz(acf)` is obtained as a biproduct of the Generalized Schur algorithm.
///
/// @return Scalar containing the value of `log(det(Tz))`.
inline double Toeplitz::log_det() {
  if (!has_solve_) solve_setup();
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    return logdet_;
  } else {
    // N = 1 case.
    return log(acf_[0]);
  }
}

/// @return The value of `trace(Tz^{-1})`, where `Tz = Toeplitz(acf)`. 
inline double Toeplitz::trace_inv() {
  if(!has_trace_) {
    if (!has_solve_) solve_setup();
    traceinv_ = 0.0;
    for (int ii = 0; ii < N_; ii++) {
      traceinv_ += (N_ - 2 * ii) * delta_[ii] * delta_[ii];
    }
    traceinv_ /= delta_[0];
    has_trace_ = true;
  }
  return traceinv_;
}

/// Computes `y = Tz * x`, where `Tz = Toeplitz(acf1)`.  This is a convenience function which does not destroy the `Toeplitz` intermediate calculations by calling `set_acf(acf1)`.
///
/// @param[out] y Vector of length `N` containing the output `y = Tz * x`.
/// @param[in] x Input vector of length `N`.
/// @param[in] acf1 Vector of length `N` specifying the first row/column of the Toeplitz matrix `Tz = Toeplitz(acf1)`.
inline void Toeplitz::prod(double* y, const double* x, const double* acf1) {
  // Pointers to temporary storage: x_, x_fft_, z_, z_fft_, y_.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* z_ = vec2_;
  dcomplex* z_fft_ = vec2_fft_;
  double* y_ = vec3_;
  // dcomplex* y_fft_ = vec2_fft_;
  // circulant embedding of acf1
  std::copy(acf1, acf1 + N_, z_);
  z_[N_] = 0;
  std::reverse_copy(acf1 + 1, acf1 + N_, z_ + N_ + 1);
  // std::copy(acf1 + 1, acf1 + N_, z_ + N_ + 1);
  // std::reverse(z_ + N_ + 1, z_ + 2 * N_);
  rfft_->fft(z_fft_, z_);
  // y_ = ifft(fft(z_) * fft(x_))[1:N_]
  std::copy(x, x + N_, x_);
  zero_fft(x_fft_, x_);
  conv_fft(y_, z_fft_, x_fft_);
  std::copy(y_, y_ + N_, y);
  return;
}

/// Computes the Toeplitz matrix-vector product `y = Tz * x`, where `Tz` is an non-symmetric Toeplitz matrix with first row `row1` and first column `col1`.  This method does not call `set_acf()`, and therefore does not reset the intermediate Toeplitz calculations.
///
/// @param[out] y Vector of length `N` containing the output `y = Tz * x`.
/// @param[in] x Input vector of length `N`.
/// @param[in] col1 Vector of length `N` specifying the first column of `Tz`.
/// @param[in] row1 Vector of length `N` specifying the first row of `Tz`.
inline void Toeplitz::prod(double* y, const double* x,
			      const double* col1, const double* row1) {
  // Pointers to temporary storage: x_, x_fft_, z_, z_fft_, y_.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* z_ = vec2_;
  dcomplex* z_fft_ = vec2_fft_;
  double* y_ = vec3_;
  // z_ = [col1, 0, rev(row1[-1])]
  std::copy(col1, col1 + N_, z_);
  z_[N_] = 0;
  std::reverse_copy(row1 + 1, row1 + N_, z_ + N_ + 1);
  // std::reverse(z_ + N_ + 1, z_ + 2 * N_);
  rfft_->fft(z_fft_, z_);
  // y_ = ifft(fft(z_) * fft(x_))[1:N_]
  std::copy(x, x + N_, x_);
  zero_fft(x_fft_, x_);
  conv_fft(y_, z_fft_, x_fft_);
  std::copy(y_, y_ + N_, y);
  return;
}

/// Computes `trace( Tz^{-1} * Tz0 )`, where `Tz = Toeplitz(acf)` and `Tz0 = Toeplitz(acf0)`.  This trace-product appears in the computation of
/// ```
///  d/dx log(det(Tz)) = trace( Tz^{-1} * Toeplitz(d/dx Tz),
/// ```
/// i.e., where `Tz = Toeplitz(acf(x))` is a function of `x`.
///
/// @param[in] acf0 Vector of length `N` giving the first row/column of the Toeplitz matrix `Tz0 = Toeplitz(acf0)`.
/// @return The Toeplitz trace-product `trace( Tz^{-1} * Tz0 )`.
inline double Toeplitz::trace_grad(const double* acf0) {
  // Pointers to temporary storage: U1_, U1_fft_, U2_, U2_fft_, y_.
  double* U1_ = vec1_;
  dcomplex* U1_fft_ = vec1_fft_;
  double* U2_ = vec2_;
  dcomplex* U2_fft_ = vec2_fft_;
  double* y_ = vec3_;
  // dcomplex* y_fft_ = vec3_fft_;
  double trace;
  double acf00 = acf0[0];
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    if (!has_solve_) solve_setup();
    // check first term to avoid singularity
    bool sng = fabs(acf0[0]) < 0.0001;
    if(sng) acf00 += 1.0;
    std::copy(acf0, acf0 + N_, U1_);
    if(sng) U1_[0] += 1.0;
    zero_fft(U1_fft_, U1_);
    std::fill(U2_, U2_ + 2 * N_, 0.0);
    std::copy(acf0 + 1, acf0 + N_, U2_ + 1);
    rfft_->fft(U2_fft_, U2_);  // U2
    // tr{U1'L1L1'U1}
    conv_fft(y_, L1_fft_, U1_fft_);
    trace = trace_LU(y_, y_);
    // tr{U1'L2L2'U1}
    conv_fft(y_, L2_fft_, U1_fft_);
    trace -= trace_LU(y_, y_);
    // tr{U2'L1L1'U2}
    conv_fft(y_, L1_fft_, U2_fft_);
    trace -= trace_LU(y_, y_);
    // tr{U2'L2L2'U2}
    conv_fft(y_, L2_fft_, U2_fft_);
    trace += trace_LU(y_, y_);
    // trace
    trace /= delta_[0];
    trace /= acf00;
    if (sng) trace -= trace_inv(); // singularity correction
  }
  else {
    trace = acf0[0] / acf_[0]; // N = 1 case.
  }
  return trace;
}

/// Computes `trace((Tz^{-1} * Tz1) * (Tz^{-1} * Tz2))`, where `Tz = Toeplitz(acf)`, `Tz1 = Toeplitz(acf1)`, and `Tz2 = Toeplitz(acf2)`.  This trace-product appears in the computation of
/// ```
/// d^2/dxdy log(det(Tz)) = trace(Tz^{-1} d^2/dxdy Tz) - trace((Tz^{-1} d/dx Tz) * (Tz^{-1} d/dy Tz)),
/// ```
/// i.e., where `Tz = Toeplitz(acf(x,y))` is a function of `x` and `y`.
///
/// @param[in] acf1 Vector of length `N` giving the first row/column of the  Toeplitz matrix `Tz1 = Toeplitz(acf1)`.
/// @param[in] acf2 Vector of length `N` giving the first row/column of the  Toeplitz matrix `Tz2 = Toeplitz(acf2)`.
/// @return The Toeplitz trace-product `trace((Tz^{-1} * Tz1) * (Tz^{-1} * Tz2))`.
inline double Toeplitz::trace_hess(const double* acf1, const double* acf2) {
  // Pointers to temporaries: phi_, x_, x_fft_, y_, y_fft_, z_, z_fft_, U1_, U2_
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  dcomplex* y_fft_ = vec2_fft_;
  double* z_ = vec3_;
  dcomplex* z_fft_ = vec3_fft_;
  double* U1_ = vec4_;
  // dcomplex* U1_fft_ = vec4_fft_;
  double* U2_ = vec5_;
  // dcomplex* U2_fft_ = vec5_fft_;
  double* phi_ = vec6_;
  double trace, kappa1, kappa2;
  double acf20 = acf2[0];
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    if (!has_solve_) solve_setup();
    // check first term to avoid singularity
    bool sng = fabs(acf2[0]) < 0.0001;
    if(sng) acf20 += 1.0;		
    // Store the negative derivative of delta in vector phi_, where phi_ = solve(acf_) * toep(acf1) * delta
    prod(phi_, delta_, acf1);
    solve(phi_, phi_);
    trace = trace_grad(acf2);
    if(sng) trace += trace_inv();
    trace *= -phi_[0];
    // kappa1
    // U1_ = phi_ \conv acf2
    std::copy(phi_, phi_ + N_, x_);
    zero_fft(x_fft_, x_);
    std::copy(acf2, acf2 + N_, z_);
    if(sng) z_[0] += 1.0;
    zero_fft(z_fft_, z_);
    conv_fft(U1_, x_fft_, z_fft_);
    // U2_ = delta_ \conv acf2
    std::copy(delta_, delta_ + N_, y_);
    zero_fft(y_fft_, y_);
    conv_fft(U2_, y_fft_, z_fft_);
    kappa1 = trace_LU(U1_, U2_) / acf20;
    // U1_ = phi_ \conv zacf2, where zacf2 = acf2 with zero in first entry
    z_[0] = 0.0;
    rfft_->fft(z_fft_, z_);
    conv_fft(U1_, x_fft_, z_fft_);
    // U2_ = delta_ \conv zacf2
    conv_fft(U2_, y_fft_, z_fft_);
    kappa1 -= trace_LU(U1_, U2_) / acf20;
    // kappa2
    // U1_ = zrphi \conv acf2, where zrphi = phi with zero in first entry, then reversed
    x_[0] = 0.0;
    std::reverse(x_ + 1, x_ + N_);
    rfft_->fft(x_fft_, x_); 
    z_[0] = acf20;
    rfft_->fft(z_fft_, z_);
    conv_fft(U1_, x_fft_, z_fft_);
    // U2_ = zrdelta \conv acf2
    y_[0] = 0.0;
    std::reverse(y_ + 1, y_ + N_);
    rfft_->fft(y_fft_, y_); 
    conv_fft(U2_, y_fft_, z_fft_);
    kappa2 = trace_LU(U1_, U2_) / acf20;
    // U1_ = zrphi \conv zacf2
    z_[0] = 0.0;
    rfft_->fft(z_fft_, z_); 
    conv_fft(U1_, x_fft_, z_fft_);
    // U2_ = zrdelta \conv zacf2
    conv_fft(U2_, y_fft_, z_fft_);
    kappa2 -= trace_LU(U1_, U2_) / acf20;
    // finalize trace calculation
    trace += 2.0 * (kappa1 - kappa2);
    trace /= delta_[0];
    if (sng) {
      // singularity correction
      double t0 = 0.0;
      for (int ii = 0; ii < N_; ii++) {
	t0 += (N_ - 2 * ii) * delta_[ii] * phi_[ii];
      }
      trace -= 2.0 * t0 / delta_[0];
      trace += trace_inv() * phi_[0] / delta_[0];
    }		
  }
  else {
    trace = acf1[0] * acf2[0] / acf_[0] / acf_[0]; // N = 1 case.
  }
  return trace;
}

#endif

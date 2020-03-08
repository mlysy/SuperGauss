/// @file Toeplitz.h

///////////////////////////////////////////////
// Toeplitz matrix class
///////////////////////////////////////////////

#ifndef Toeplitz_h
#define Toeplitz_h 1

#include "GSchur.h"

/// Trace of produce of lower and upper triangular Toeplitz matrices.
///
/// Calculates `trace(L * U)`, where `L` and `U` are lower/upper triangular Toeplitz matrices. It is an `O(n)` algorithm.
///
/// @param[in] L First column of `L`.
/// @param[in] U First row of `U`.
/// @param[in] n Size of inputs.
inline double trace_LU(double* L, double* U, int n) {
  double trace = 0;
  for (int ii = 0; ii < n; ++ii) {
    trace += (n - ii) * L[ii] * U[ii];
  }
  return trace;
}

/// Class for Toeplitz-related computation.
///
/// Generate the memory for Toeplitz computation including multiplication,
/// inversion, determinant, first and second order derivative of the
/// log-determinant.
class Toeplitz {
protected:
  typedef std::complex<double> dcomplex;
  // internal storage.  don't overwrite locally.
  int N_;        ///< Size of Toeplitz matrix.
  int N2_;       ///< Size of FFT product.
  int N3_;       ///< Size of other FFT product.
  double* acf_;    ///< First column of the Toeplitz matrix.
  double* tzcirc_; ///< Storage for Toeplitz circulant embedding.
  dcomplex* tzcirc_fft_; ///< FFT of Toeplitz circulant embedding.
  double* delta_; ///< Storage for first column of the inverse-Toeplitz matrix.
  double logdet_; ///< Storage for log-determinant of Toeplitz matrix.
  double traceinv_;///< Trace of inverse-Toeplitz.
  dcomplex* vec_fft_; ///< FFT for convolutions.
  double* vec_zero_; ///< Storage for zero-padding
  GSchurN* gs_;     ///< Object to compute GSchur algorithm.
  VectorFFT* vfft_;   ///< Object for fft computations.
  bool has_acf_;    ///< Boolean flag indicating whether input argument `acf_` has been modified.
  bool has_mult_;   ///< Boolean flag indicating whether multiplication-related FFT has been done.
  bool has_solve_;  ///< Boolean flag indicating whether inversion-related FFT has been done.
  bool has_trace_;  ///< Boolean flag indicating whether trace-of-inverse calculation has been done.
  dcomplex* L1_fft_; ///< FFT for first lower-triangular Toeplitz matrix in `solve`.
  dcomplex* tL1_fft_; ///< FFT for transpose of first lower-triangular Toeplitz matrix in `solve`.
  dcomplex* L2_fft_; ///< FFT for second lower-triangular Toeplitz matrix in `solve`.
  dcomplex* tL2_fft_; ///< FFT for transpose of second lower-triangular Toeplitz matrix in `solve`.
  // temporary storage.  ok to overwrite locally.
  double *vec1_, *vec2_, *vec3_, *vec4_, *vec5_, *vec6_;
  dcomplex *vec1_fft_, *vec2_fft_, *vec3_fft_, *vec4_fft_, *vec5_fft_;
  // double* x_;
  // double* y_;      // Storage for `mult/solve` output.
  // double* z_;      // Storage for `mult/solve/trace` temporary.
  // double* phi_;    // Storage for `trace` temporary.
  // double* U1_;     // Storage for first upper triangular toeplitz matrix.
  // double* U2_;     // Storage for second upper triangular toeplitz matrix.
  // dcomplex* x_fft_;
  // dcomplex* y_fft_;
  // dcomplex* z_fft_;
  // dcomplex* U1_fft_;
  // dcomplex* U2_fft_;  
  /// Prepare for multiplication.
  void mult_setup();    
  /// Prepare for solving linear systems.
  void solve_setup();
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
  void setAcf(const double* acf);
  /// Get the acf of the Toepliz matrix.
  void getAcf(double* acf);  
  /// Size of the Toeplitz matrix.
  int size(); 
  /// Check whether the acf has been set.
  bool hasAcf();  
  /// Toeplitz matrix-vector multiplication.
  void mult(double* y, const double* x);
  /// External symmetric Toeplitz matrix-vector multiplication.
  void product(double* y, const double* x, const double* acf1);
  /// External asymmetric Toeplitz matrix-vector multiplication.
  void product(double* y, const double* x, const double* acf1, const double* acf2);  
  /// Solve Toeplitz matrix-vector system of equations.
  void solve(double* y, const double* x);
  /// Log-determinant of Toeplitz matrix.
  double logDet();
  /// Trace of inverse-Toeplitz matrix.
  double trace_inv();
  /// Gradient-specialized trace-product.
  double trace_deriv(const double* acf1);
  /// Hessian-specialized trace-product.
  double trace_hess(const double* acf1, const double* acf2);
};

/// @param[in] N Size of Toeplitz matrix.
/// @param[in] bmod Size of binary modulus for GSchur calculation.
inline Toeplitz::Toeplitz(int N, int bmod = 5) {
  N_ = N;
  N2_ = 2 * (N_ / 2 + 1);
  N3_ = N2_; // N_ + 1;
  acf_ = new double[N_];
  has_acf_ = false;
  has_mult_ = false;
  has_solve_ = false;
  has_trace_ = false;
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1
    gs_ = new GSchurN(N_, bmod);
    vfft_ = new VectorFFT(2 * N_);
    tzcirc_ = new double[2 * N_];
    tzcirc_fft_ = new dcomplex[2 * N_];
    vec_fft_ = new dcomplex[2 * N_];
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
    vec4_fft_ = new dcomplex[2 * N_];
    vec5_ = new double[2 * N_];
    vec5_fft_ = new dcomplex[2 * N_];
    vec6_ = new double[2 * N_];
  }
}

inline Toeplitz::~Toeplitz() {
  delete[] acf_;
  // GSchur algorithm only supports N > 1 case.
  if (N_ > 1) {
    delete gs_;
    delete vfft_;
    delete[] tzcirc_;
    delete[] tzcirc_fft_;
    delete[] vec_fft_;
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
    delete[] vec4_fft_;
    delete[] vec5_;
    delete[] vec5_fft_;
    delete[] vec6_;    
    // delete[] x_;
    // delete[] x_fft_;
    // delete[] y_;
    // delete[] y_fft_;
    // delete[] z_;
    // delete[] z_fft_;
    // delete[] phi_;    
    // delete[] U1_;
    // delete[] U2_;
    // delete[] U1_fft_;
    // delete[] U2_fft_;
  }
}


/// @param[in] acf Real vector of ACF.
inline void Toeplitz::setAcf(const double* acf) {
  std::copy(acf, acf + N_, acf_);
  has_acf_ = true;
  has_mult_ = false;
  has_solve_ = false;
  has_trace_ = false;
  return;
}

/// @param[out] acf Real vector of ACF.
inline void Toeplitz::getAcf(double* acf) {
  std::copy(acf_, acf_ + N_, acf);
  return;
}

inline bool Toeplitz::hasAcf() { 
  return has_acf_; 
}

inline int Toeplitz::size() { 
  return N_; 
}

/// Preparation of computing \f$T \times x\f$, i.e. generating
/// FFT([acf_, 0, rev(acf_[-1])]) and change the flag `has_mult_`.
/// FFT result is stored in tzcirc_fft_.
inline void Toeplitz::mult_setup() {
  has_mult_ = true;
  if (N_ > 1) {
    std::copy(acf_, acf_ + N_, tzcirc_);
    tzcirc_[N_] = 0;
    std::copy(acf_ + 1, acf_ + N_, tzcirc_ + N_ + 1);
    std::reverse(tzcirc_ + N_ + 1, tzcirc_ + 2 * N_);
    vfft_->fft(tzcirc_fft_, tzcirc_);
  }
  return;
}

/// Calculates the convolution between length-`N_` vectors `x` and `y` based on their (zero-padded) FFT transformations `x_fft` and `y_fft`.  This is an elementwise multiplication of `x_fft` and `y_fft` followed by an iFFT back to the real domain.
///
/// @param[out] z Convolution output.
/// @param[in] x_fft FFT of first input vector.
/// @param[in] y_fft FFT of second input vector.
inline void Toeplitz::conv_fft(double* z, const dcomplex* x_fft, const dcomplex* y_fft) {
  complex_mult(vec_fft_, x_fft, y_fft, N2_);
  vfft_->ifft(z, vec_fft_);
  return;
}

/// Calculates the FFT of the length-`2N_` vector `[x, 0, ..., 0]`.
///
/// @param[out] y_fft FFT of zero-padded vector `x`.
/// @param[in] x Input vector.
///
/// @warning Last `N_` elements of input vector `x` are _modified_.
inline void Toeplitz::zero_fft(dcomplex* y_fft, double* x) {
  std::fill(x + N_, x + 2 * N_, 0.0);
  vfft_->fft(y_fft, x);
}

/// Product \f$ T \times x\f$ is computed efficiently by extending Toeplitz
/// matrix \f$T\f$ into a circulant matrix, whose multiplication with vector can
/// be done with Fast Fourier Transformation (FFT) in \f$O(n \log n) \f$ steps.
///
/// @param[out] y Real vector of result.
/// @param[in] x Real vector for product.
inline void Toeplitz::mult(double* y, const double* x) {
  // Pointers to temporary storage: x_, x_fft_, y_.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  // dcomplex* y_fft_ = vec2_fft_;
  // tzcirc_ = [acf_, 0, rev(acf_[-1])]
  if(!has_mult_) mult_setup();
  // x_ = [x, 0]
  std::copy(x, x + N_, x_);
  // std::fill(x_ + N_, x_ + 2 * N_, 0.0);
  // vfft_->fft(x_fft_, x_);
  zero_fft(x_fft_, x_);
  // y_ = ifft(fft(tzcirc_) * fft(x_))[1:N_]
  // complex_mult(y_fft_, tzcirc_fft_, x_fft_, N2_);
  // vfft_->ifft(y_, y_fft_);
  conv_fft(y_, tzcirc_fft_, x_fft_);
  std::copy(y_, y_ + N_, y);
  return;
}



/// Preparation of computing \f$T^{-1} \times x\f$, including GSchur algorithm
/// implementation and generating FFT for Gohberg-Semencul formula.
inline void Toeplitz::solve_setup() {
  // Pointers to temporary storage: `z_`.
  double* z_ = vec1_;
  has_solve_ = true;
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    gs_->compute(delta_, logdet_, acf_);
    /// tL1_fft_ stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_1'
    z_[0] = delta_[0];
    std::fill(z_ + 1, z_ + N_ + 1, 0);
    // std::copy(delta_ + 1, delta_ + N_, z_ + N_ + 1);
    // std::reverse(z_ + N_ + 1, z_ + 2 * N_);
    std::reverse_copy(delta_ + 1, delta_ + N_, z_ + N_ + 1);
    vfft_->fft(tL1_fft_, z_);
    /// L1_fft_ stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_1
    std::copy(delta_, delta_ + N_, z_);
    // std::fill(z_ + N_, z_ + 2 * N_, 0);
    // vfft_->fft(L1_fft_, z_);
    zero_fft(L1_fft_, z_);
    /// tL2_fft_ stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_2'
    std::fill(z_, z_ + N_ + 1, 0);
    std::copy(delta_ + 1, delta_ + N_, z_ + N_ + 1);
    vfft_->fft(tL2_fft_, z_);	
    /// L2_fft_ stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_2
    std::fill(z_, z_ + 2 * N_, 0);
    // std::copy(delta_ + 1, delta_ + N_, z_ + 1);
    // std::reverse(z_ + 1, z_ + N_);
    std::reverse_copy(delta_ + 1, delta_ + N_, z_ + 1);
    vfft_->fft(L2_fft_, z_);
  }
  return;
}

/// Multiplication between inverse Toeplitz matrix and vector: \f$T^{-1} \times
/// x\f$ involves the Gohberg-Semencul formula which decomposes inverse Toeplitz
/// into sum of products of lower and upper triangular Toeplitz matrices: \f$
/// T^{-1} = \frac{1}{\rho_0} [L_1L_1' - L_2 L_2'] \f$, where \f$L_1, L_2\f$ are
/// lower-trangular Toeplitz matrices and \f$\rho_0\f$ is the first element of
/// GSchur Algorithm output: length-N vector \f$\rho\f$.
///
/// @param[out] y Real vector of result.
/// @param[in] x Real vector for product.
inline void Toeplitz::solve(double* y, const double* x) {
  // Pointers to temporary storage: `x_`, `x_fft_`, `y_`, `y_fft_`, `z_`, `z_fft_`.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  dcomplex* y_fft_ = vec2_fft_;
  double* z_ = vec3_;
  dcomplex* z_fft_ = vec3_fft_;
  if (!has_solve_) solve_setup();
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(x, x + N_, x_);
    // std::fill(x_ + N_, x_ + 2 * N_, 0);
    // vfft_->fft(x_fft_, x_);
    zero_fft(x_fft_, x_);
    // y_ = L1' * x_
    // complex_mult(y_fft_, tL1_fft_, x_fft_, N3_);
    // vfft_->ifft(y_, y_fft_);
    conv_fft(y_, tL1_fft_, x_fft_);
    // z_ = L1 * y_ = L1 * L1' * x_
    // std::fill(y_ + N_, y_ + 2 * N_, 0);
    // vfft_->fft(y_fft_, y_);
    zero_fft(y_fft_, y_);
    // complex_mult(z_fft_, L1_fft_, y_fft_, N3_);
    // vfft_->ifft(z_, z_fft_);
    conv_fft(z_, L1_fft_, y_fft_);
    // y_ = L2' * x_
    // complex_mult(y_fft_, tL2_fft_, x_fft_, N3_);
    // vfft_->ifft(y_, y_fft_);
    conv_fft(y_, tL2_fft_, x_fft_);
    // x_ = L2 * y_ = L2 * L2' * x_ (temporarily stored in x_)
    // std::fill(y_ + N_, y_ + 2 * N_, 0);
    // vfft_->fft(y_fft_, y_);
    zero_fft(y_fft_, y_);
    // complex_mult(x_fft_, L2_fft_, y_fft_, N3_);
    // vfft_->ifft(x_, x_fft_);
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

/// Log-determinant of Toeplitz matrix.
inline double Toeplitz::logDet() {
  if (!has_solve_) solve_setup();
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    return logdet_;
  } else {
    // N = 1 case.
    return log(acf_[0]);
  }
}

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

/// Computing the multiplication between a symmetric Toeplitz matrix and a vector.
/// This is a convenient function that takes outside arguments including acf1, acf2 and x and has nothing to do with the built-in argument acf_.
/// We kind of borrow the memory of Toeplitz class for such computation.
/// @param[out] y Real vector of result.
/// @param[in] x Real vector for product.
/// @param[in] acf1 Real vector of the first column of the asymmetric Toeplitz matrix.
inline void Toeplitz::product(double* y, const double* x, const double* acf1) {
  // Pointers to temporary storage: `x_`, `x_fft_`, `y_`, `y_fft_`, `z_`, `z_fft_`.
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  dcomplex* y_fft_ = vec2_fft_;
  double* z_ = vec3_;
  dcomplex* z_fft_ = vec3_fft_;
  // circulant embedding of acf1
  std::copy(acf1, acf1 + N_, z_);
  z_[N_] = 0;
  std::copy(acf1 + 1, acf1 + N_, z_ + N_ + 1);
  std::reverse(z_ + N_ + 1, z_ + 2 * N_);
  vfft_->fft(z_fft_, z_);
  // x_ = [x, 0]
  std::copy(x, x + N_, x_);
  // std::fill(x_ + N_, x_ + 2 * N_, 0);
  // vfft_->fft(x_fft_, x_);
  zero_fft(x_fft_, x_);
  // y_ = ifft(fft(z_) * fft(x_))[1:N_]
  // complex_mult(y_fft_, z_fft_, x_fft_, N3_);
  // vfft_->ifft(y_, y_fft_);
  conv_fft(y_, z_fft_, x_fft_);
  std::copy(y_, y_ + N_, y);
  return;
}

/// Computing the multiplication between a non-symmetric Toeplitz matrix and a vector.
/// This is a convenient function that takes outside arguments including acf1, acf2 and x and has nothing to do with the built-in argument acf_.
/// We kind of borrow the memory of Toeplitz class for such computation.
/// @param[out] y Real vector of result.
/// @param[in] x Real vector for product.
/// @param[in] acf1 Real vector of the first column of the asymmetric Toeplitz matrix.
/// @param[in] acf2 Real vector of the first row of the asymmetric Toeplitz matrix.
inline void Toeplitz::product(double* y, const double* x, const double* acf1, const double* acf2) {
  // Pointers to temporary storage: `x_`, `x_fft_`, `y_`, `y_fft_`, `z_`, `z_fft_`. 
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  dcomplex* y_fft_ = vec2_fft_;
  double* z_ = vec3_;
  dcomplex* z_fft_ = vec3_fft_;
  // z_ = [acf1, 0, rev(acf2[-1])]
  std::copy(acf1, acf1 + N_, z_);
  z_[N_] = 0;
  std::copy(acf2 + 1, acf2 + N_, z_ + N_ + 1);
  std::reverse(z_ + N_ + 1, z_ + 2 * N_);
  vfft_->fft(z_fft_, z_);
  // x_ = [x, 0]
  std::copy(x, x + N_, x_);
  // std::fill(x_ + N_, x_ + 2 * N_, 0);
  // vfft_->fft(x_fft_, x_);
  zero_fft(x_fft_, x_);
  // y_ = ifft(fft(z_) * fft(x_))[1:N_]
  // complex_mult(y_fft_, z_fft_, x_fft_, N2_);
  // vfft_->ifft(y_, y_fft_);
  conv_fft(y_, z_fft_, x_fft_);
  std::copy(y_, y_ + N_, y);
  return;
}


inline double Toeplitz::trace_deriv(const double* acf0) {
  // Pointers to temporary storage: `U1_`, `U1_fft_`, `U2_`, `U2_fft_`, `y_`, `y_fft_`.
  double* U1_ = vec1_;
  dcomplex* U1_fft_ = vec1_fft_;
  double* U2_ = vec2_;
  dcomplex* U2_fft_ = vec2_fft_;
  double* y_ = vec3_;
  dcomplex* y_fft_ = vec3_fft_;
  double trace;
  double acf00 = acf0[0];
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    if (!has_solve_) solve_setup();
    // check first term to avoid singularity
    bool sng = fabs(acf0[0] < 0.0001);
    if(sng) acf00 += 1.0;
    // bool sng = false;
    // if (fabs(acf0[0]) < 0.0001) {
    //   sng = true;
    //   acf0[0] += 1;
    // }
    std::copy(acf0, acf0 + N_, U1_);
    if(sng) U1_[0] += 1.0;
    // std::fill(U1_ + N_, U1_ + 2 * N_, 0);
    // vfft_->fft(U1_fft_, U1_); // U1
    zero_fft(U1_fft_, U1_);
    std::fill(U2_, U2_ + 2 * N_, 0);
    std::copy(acf0 + 1, acf0 + N_, U2_ + 1);
    vfft_->fft(U2_fft_, U2_);  // U2
    // tr{U1'L1L1'U1}
    // complex_mult(y_fft_, L1_fft_, U1_fft_, N3_);
    // vfft_->ifft(y_, y_fft_);
    conv_fft(y_, L1_fft_, U1_fft_);
    trace = trace_LU(y_, y_, N_);
    // tr{U1'L2L2'U1}
    // complex_mult(y_fft_, L2_fft_, U1_fft_, N2_);
    // vfft_->ifft(y_, y_fft_);
    conv_fft(y_, L2_fft_, U1_fft_);
    trace -= trace_LU(y_, y_, N_);
    // tr{U2'L1L1'U2}
    // complex_mult(y_fft_, L1_fft_, U2_fft_, N2_);
    // vfft_->ifft(y_, y_fft_);
    conv_fft(y_, L1_fft_, U2_fft_);
    trace -= trace_LU(y_, y_, N_);
    // tr{U2'L2L2'U2}
    // complex_mult(y_fft_, L2_fft_, U2_fft_, N2_);
    // vfft_->ifft(y_, y_fft_);
    conv_fft(y_, L2_fft_, U2_fft_);
    trace += trace_LU(y_, y_, N_);
    // trace
    trace /= delta_[0];
    trace /= acf00;
    if (sng) {
      trace -= trace_inv();
      // double t0 = 0.0;
      // for (int ii = 0; ii < N_; ii++) {
      // 	t0 += (N_ - 2 * ii) * delta_[ii] * delta_[ii];
      // }
      // trace -= t0 / delta_[0];
      // acf0[0] -= 1;

    }
  }
  else {
    // N = 1 case.
    trace = acf0[0] / acf_[0];
  }
  return trace;
}

/// Function for another trace computation
inline double Toeplitz::trace_hess(const double* acf1, const double* acf2) {
  // Pointers to temporaries: phi_, x_, x_fft_, y_, y_fft_, z_, z_fft_, U1_, U1_fft, U2_, U2_fft_
  double* x_ = vec1_;
  dcomplex* x_fft_ = vec1_fft_;
  double* y_ = vec2_;
  dcomplex* y_fft_ = vec2_fft_;
  double* z_ = vec3_;
  dcomplex* z_fft_ = vec3_fft_;
  double* U1_ = vec4_;
  dcomplex* U1_fft_ = vec4_fft_;
  double* U2_ = vec5_;
  dcomplex* U2_fft_ = vec5_fft_;
  double* phi_ = vec6_;
  double trace, kappa1, kappa2;
  double acf20 = acf2[0];
  if (N_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    if (!has_solve_) solve_setup();
    // bool sng = false;
    // if (fabs(acf2[0]) < 0.0001) {
    //   sng = true;
    //   acf2[0] += 1;
    // }
    // check first term to avoid singularity
    bool sng = fabs(acf2[0] < 0.0001);
    if(sng) acf20 += 1.0;		
    // Store the negative derivative of delta in vector phi_, where phi_ = solve(acf_) * toep(acf1) * delta
    product(phi_, delta_, acf1);
    solve(phi_, phi_);
    trace = trace_deriv(acf2);
    if(sng) trace += trace_inv();
    trace *= -phi_[0];
    // kappa1
    std::copy(phi_, phi_ + N_, x_);
    // std::fill(x_ + N_, x_ + 2 * N_, 0);
    // vfft_->fft(x_fft_, x_); // x_fft_ = fft(phi_zero_)
    zero_fft(x_fft_, x_);
    std::copy(acf2, acf2 + N_, z_);
    if(sng) z_[0] += 1.0;
    // std::fill(z_ + N_, z_ + 2 * N_, 0.0);
    // vfft_->fft(z_fft_, z_); // z_fft_ = fft(acf2_zero_)
    zero_fft(z_fft_, z_);
    // complex_mult(U1_fft_, x_fft_, z_fft_, N3_);
    // vfft_->ifft(U1_, U1_fft_); // U1_ = phi_ \conv acf2
    conv_fft(U1_, x_fft_, z_fft_);
    std::copy(delta_, delta_ + N_, y_);
    // std::fill(y_ + N_, y_ + 2 * N_, 0);
    // vfft_->fft(y_fft_, y_); // y_fft_ = fft(delta_zero)
    zero_fft(y_fft_, y_);
    // complex_mult(U2_fft_, y_fft_, z_fft_, N3_);
    // vfft_->ifft(U2_, U2_fft_); // U2_ = delta_ \conv acf2
    conv_fft(U2_, y_fft_, z_fft_);
    kappa1 = trace_LU(U1_, U2_, N_) / acf20;
    z_[0] = 0.0; // z_ = zacf2 = acf2 with zero in first entry
    vfft_->fft(z_fft_, z_);
    // complex_mult(U1_fft_, x_fft_, z_fft_, N3_);
    // vfft_->ifft(U1_, U1_fft_); // U1_ = phi_ \conv zacf2
    conv_fft(U1_, x_fft_, z_fft_);
    // complex_mult(U2_fft_, y_fft_, z_fft_, N3_);
    // vfft_->ifft(U2_, U2_fft_); // U2_ = delta_ \conv zacf2
    conv_fft(U2_, y_fft_, z_fft_);
    kappa1 -= trace_LU(U1_, U2_, N_) / acf20;
    // kappa2
    x_[0] = 0;
    std::reverse(x_ + 1, x_ + N_);
    vfft_->fft(x_fft_, x_); // x_fft_ = fft(zrphi), where rzphi = phi with zero in first entry, then reversed
    z_[0] = acf20;
    vfft_->fft(z_fft_, z_); // z_fft_ = fft(acf2)
    // complex_mult(U1_fft_, x_fft_, z_fft_, N3_);
    // vfft_->ifft(U1_, U1_fft_); // U1_ = zrphi \conv acf2
    conv_fft(U1_, x_fft_, z_fft_);
    y_[0] = 0;
    std::reverse(y_ + 1, y_ + N_);
    vfft_->fft(y_fft_, y_); // y_fft_ = fft(zrdelta)
    // complex_mult(U2_fft_, y_fft_, z_fft_, N3_);
    // vfft_->ifft(U2_, U2_fft_); // U2_ = zrdelta \conv acf2
    conv_fft(U2_, y_fft_, z_fft_);
    kappa2 = trace_LU(U1_, U2_, N_) / acf20;
    z_[0] = 0;
    vfft_->fft(z_fft_, z_); // z_fft_ = fft(zacf2)
    // complex_mult(U1_fft_, x_fft_, z_fft_, N3_);
    // vfft_->ifft(U1_, U1_fft_); // U1_ = zrphi \conv zacf2
    conv_fft(U1_, x_fft_, z_fft_);
    // complex_mult(U2_fft_, y_fft_, z_fft_, N3_);
    // vfft_->ifft(U2_, U2_fft_); // U2_ = zrdelta \conv zacf2
    conv_fft(U2_, y_fft_, z_fft_);
    kappa2 -= trace_LU(U1_, U2_, N_) / acf20;
    // finalize trace calculation
    trace += 2 * (kappa1 - kappa2);
    trace /= delta_[0];
    if (sng) {
      double t0 = 0.0;
      for (int ii = 0; ii < N_; ii++) {
	t0 += (N_ - 2 * ii) * delta_[ii] * phi_[ii];
      }
      trace -= 2 * t0 / delta_[0];
      trace += trace_inv() * phi_[0] / delta_[0];
      // t0 = 0.0;
      // for (int ii = 0; ii < N_; ii++) {
      // 	t0 += (N_ - 2 * ii) * delta_[ii] * delta_[ii];
      // }
      // trace += t0 * phi_[0] / delta_[0] / delta_[0];
      // acf2[0] -= 1;
    }		
  }
  else {
    // N = 1 case.
    trace = acf1[0] * acf2[0] / acf_[0] / acf_[0];
  }
  return trace;
}


#endif

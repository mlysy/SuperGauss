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
  int ntz_;        ///< Size of Toeplitz matrix.
  double* acf_;    ///< First column of the Toeplitz matrix.
  double* tzcirc_; ///< Storage for Toeplitz circulant embedding.
  dcomplex* tzcirc_fft_; ///< FFT of Toeplitz circulant embedding.
  GSchurN* gs_;     ///< Object to compute GSchur algorithm.
  VectorFFT* vfft_;   ///< Object for fft computations.
  double traceinv_;///< Trace of inverse-Toeplitz.
  bool has_acf_;    ///< Boolean flag indicating whether input argument `acf_` has been modified.
  bool has_mult_;   ///< Boolean flag indicating whether multiplication-related FFT has been done.
  bool has_solve_;  ///< Boolean flag indicating whether inversion-related FFT has been done.
  bool has_trace_;  ///< Boolean flag indicating whether trace-of-inverse calculation has been done.
  dcomplex* L1_fft_; ///< FFT for first lower-triangular Toeplitz matrix in `solve`.
  dcomplex* tL1_fft_; ///< FFT for transpose of first lower-triangular Toeplitz matrix in `solve`.
  dcomplex* L2_fft_; ///< FFT for second lower-triangular Toeplitz matrix in `solve`.
  dcomplex* tL2_fft_; ///< FFT for transpose of second lower-triangular Toeplitz matrix in `solve`.
  // temporary storage.  ok to overwrite locally.
  double* x_;      // Storage for `mult/solve` input.
  double* y_;      // Storage for `mult/solve` output.
  double* z_;      // Storage for `mult/solve/trace` temporary.
  double* phi_;    // Storage for `trace` temporary.
  double* U1_;     // Storage for first upper triangular toeplitz matrix.
  double* U2_;     // Storage for second upper triangular toeplitz matrix.
  dcomplex* x_fft_;
  dcomplex* y_fft_;
  dcomplex* z_fft_;
  dcomplex* U1_fft_;
  dcomplex* U2_fft_;  
  /// Prepare for multiplication.
  void mult_setup();    
  /// Prepare for solving linear systems.
  void solve_setup();
public:
  /// Constructor.
  Toeplitz(int N);  
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
inline Toeplitz::Toeplitz(int N) {
  ntz_ = N;
  acf_ = new double[ntz_];
  has_acf_ = false;
  has_mult_ = false;
  has_solve_ = false;
  has_trace_ = false;
  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1
    gs_ = new GSchurN(ntz_);
    vfft_ = new VectorFFT(2 * ntz_);
    tzcirc_ = new double[2 * ntz_]; /// storage for circulant embedding
    tzcirc_fft_ = new dcomplex[2 * ntz_];
    L1_fft_ = new dcomplex[2 * ntz_];
    tL1_fft_ = new dcomplex[2 * ntz_];
    L2_fft_ = new dcomplex[2 * ntz_];
    tL2_fft_ = new dcomplex[2 * ntz_];
    x_ = new double[2 * ntz_]; /// storage for x
    x_fft_ = new dcomplex[2 * ntz_];
    y_ = new double[2 * ntz_]; /// storage for y
    y_fft_ = new dcomplex[2 * ntz_];
    z_ = new double[2 * ntz_]; /// storage for temporary vector
    z_fft_ = new dcomplex[2 * ntz_];
    phi_ = new double[ntz_];
    U1_ = new double[2 * ntz_];
    U2_ = new double[2 * ntz_];
    U1_fft_ = new dcomplex[2 * ntz_];
    U2_fft_ = new dcomplex[2 * ntz_];
  }
}

inline Toeplitz::~Toeplitz() {
  delete[] acf_;

  // GSchur algorithm only supports N > 1 case.
  if (ntz_ > 1) {
    delete gs_;
    delete vfft_;
    delete[] tzcirc_;
    delete[] tzcirc_fft_;
    delete[] L1_fft_;
    delete[] tL1_fft_;
    delete[] L2_fft_;
    delete[] tL2_fft_;
    delete[] x_;
    delete[] x_fft_;
    delete[] y_;
    delete[] y_fft_;
    delete[] z_;
    delete[] z_fft_;
    delete[] phi_;    
    delete[] U1_;
    delete[] U2_;
    delete[] U1_fft_;
    delete[] U2_fft_;
  }
}


/// @param[in] acf Real vector of ACF.
inline void Toeplitz::setAcf(const double* acf) {
  std::copy(acf, acf + ntz_, acf_);
  has_acf_ = true;
  has_mult_ = false;
  has_solve_ = false;
  has_trace_ = false;
  return;
}

/// @param[out] acf Real vector of ACF.
inline void Toeplitz::getAcf(double* acf) {
  std::copy(acf_, acf_ + ntz_, acf);
  return;
}

inline bool Toeplitz::hasAcf() { 
  return has_acf_; 
}

inline int Toeplitz::size() { 
  return ntz_; 
}

/// Preparation of computing \f$T \times x\f$, i.e. generating
/// FFT([acf_, 0, rev(acf_[-1])]) and change the flag `has_mult_`.
/// FFT result is stored in tzcirc_fft_.
inline void Toeplitz::mult_setup() {
  has_mult_ = true;
  if (ntz_ > 1) {
    std::copy(acf_, acf_ + ntz_, tzcirc_);
    tzcirc_[ntz_] = 0;
    std::copy(acf_ + 1, acf_ + ntz_, tzcirc_ + ntz_ + 1);
    std::reverse(tzcirc_ + ntz_ + 1, tzcirc_ + 2 * ntz_);
    vfft_->fft(tzcirc_fft_, tzcirc_);
  }
  return;
}

/// Product \f$ T \times x\f$ is computed efficiently by extending Toeplitz
/// matrix \f$T\f$ into a circulant matrix, whose multiplication with vector can
/// be done with Fast Fourier Transformation (FFT) in \f$O(n \log n) \f$ steps.
///
/// @param[out] y Real vector of result.
/// @param[in] x Real vector for product.
inline void Toeplitz::mult(double* y, const double* x) {
  /// tzcirc_ = [acf_, 0, rev(acf_[-1])]
  if(!has_mult_) mult_setup();

  /// x_ = [x, 0]
  std::copy(x, x + ntz_, x_);
  std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
  vfft_->fft(x_fft_, x_);

  // y_ = ifft(fft(tzcirc_) * fft(x_))[1:ntz_]
  complex_mult(y_fft_, tzcirc_fft_, x_fft_, 2 * (ntz_ / 2 + 1));
  vfft_->ifft(y_, y_fft_);
  std::copy(y_, y_ + ntz_, y);

  return;
}



/// Preparation of computing \f$T^{-1} \times x\f$, including GSchur algorithm
/// implementation and generating FFT for Gohberg-Semencul formula.
inline void Toeplitz::solve_setup() {
  has_solve_ = true;

  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    gs_->compute(acf_);

    /// tL1_fft_ stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_1'
    z_[0] = gs_->delta[0];
    std::fill(z_ + 1, z_ + ntz_ + 1, 0);
    std::copy(gs_->delta + 1, gs_->delta + ntz_, z_ + ntz_ + 1);
    std::reverse(z_ + ntz_ + 1, z_ + 2 * ntz_);
    vfft_->fft(tL1_fft_, z_);

    /// L1_fft_ stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_1
    std::copy(gs_->delta, gs_->delta + ntz_, z_);
    std::fill(z_ + ntz_, z_ + 2 * ntz_, 0);
    vfft_->fft(L1_fft_, z_); 

    /// tL2_fft_ stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_2'
    std::fill(z_, z_ + ntz_ + 1, 0);
    std::copy(gs_->delta + 1, gs_->delta + ntz_, z_ + ntz_ + 1);
    vfft_->fft(tL2_fft_, z_);
	
    /// L2_fft_ stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_2
    std::fill(z_, z_ + 2 * ntz_, 0);
    std::copy(gs_->delta + 1, gs_->delta + ntz_, z_ + 1);
    std::reverse(z_ + 1, z_ + ntz_);
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
  if (!has_solve_) {
    solve_setup();
  }

  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(x, x + ntz_, x_);
    std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
    vfft_->fft(x_fft_, x_);

    // y_ = L1' * x_
    complex_mult(y_fft_, tL1_fft_, x_fft_, ntz_ + 1);
    vfft_->ifft(y_, y_fft_);
    std::fill(y_ + ntz_, y_ + 2 * ntz_, 0);
    // z_ = L1 * y_ = L1 * L1' * x_
    vfft_->fft(y_fft_, y_);
    complex_mult(z_fft_, L1_fft_, y_fft_, ntz_ + 1);
    vfft_->ifft(z_, z_fft_);
    
    // y_ = L2' * x_
    complex_mult(y_fft_, tL2_fft_, x_fft_, ntz_ + 1);
    vfft_->ifft(y_, y_fft_);
    std::fill(y_ + ntz_, y_ + 2 * ntz_, 0);
    // x_ = L2 * y_ = L2 * L2' * x_ (temporarily stored in x_)
    vfft_->fft(y_fft_, y_);
    complex_mult(x_fft_, L2_fft_, y_fft_, ntz_ + 1);
    vfft_->ifft(x_, x_fft_);

    // y = (z_ - x_) / delta[1] = 1/delta[1] * (L1 * L1' * x_ - L2 * L2' * x_)
    for (int ii = 0; ii < ntz_; ++ii) {
      y[ii] = (z_[ii] - x_[ii]) / gs_->delta[0];
    }
  } else {
    // N = 1 case.
    y[0] = x[0] / acf_[0];
  }
  return;
}

/// Log-determinant of Toeplitz matrix.
inline double Toeplitz::logDet() {
  if (!has_solve_) {
    solve_setup();
  }
  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    return gs_->ldV;
  } else {
    // N = 1 case.
    return log(acf_[0]);
  }
}

inline double Toeplitz::trace_inv() {
  // double t0 = 0.0;
  if(!has_trace_) {
    if (!has_solve_) solve_setup();
    traceinv_ = 0.0;
    for (int ii = 0; ii < ntz_; ii++) {
      traceinv_ += (ntz_ - 2 * ii) * gs_->delta[ii] * gs_->delta[ii];
    }
    traceinv_ /= gs_->delta[0];
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
  std::copy(acf1, acf1 + ntz_, z_);
  z_[ntz_] = 0;
  std::copy(acf1 + 1, acf1 + ntz_, z_ + ntz_ + 1);
  std::reverse(z_ + ntz_ + 1, z_ + 2 * ntz_);
  vfft_->fft(z_fft_, z_);

  /// x_ = [x, 0]
  std::copy(x, x + ntz_, x_);
  std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
  vfft_->fft(x_fft_, x_);

  /// y_ = ifft(fft(z_) * fft(x_))[1:ntz_]
  complex_mult(y_fft_, z_fft_, x_fft_, ntz_ + 1);
  vfft_->ifft(y_, y_fft_);
  std::copy(y_, y_ + ntz_, y);

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
  /// z_ = [acf1, 0, rev(acf2[-1])]
  std::copy(acf1, acf1 + ntz_, z_);
  z_[ntz_] = 0;
  std::copy(acf2 + 1, acf2 + ntz_, z_ + ntz_ + 1);
  std::reverse(z_ + ntz_ + 1, z_ + 2 * ntz_);
  vfft_->fft(z_fft_, z_);

  /// x_ = [x, 0]
  std::copy(x, x + ntz_, x_);
  std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
  vfft_->fft(x_fft_, x_);

  /// y_ = ifft(fft(z_) * fft(x_))[1:ntz_]
  complex_mult(y_fft_, z_fft_, x_fft_, 2 * (ntz_ / 2 + 1));
  vfft_->ifft(y_, y_fft_);
  std::copy(y_, y_ + ntz_, y);

  return;
}


inline double Toeplitz::trace_deriv(const double* acf0) {
  double trace;
  double acf00 = acf0[0];

  if (ntz_ > 1) { // GSchur algorithm only supports N > 1 case.
    if (!has_solve_) {
      solve_setup();
    }

    // check first term to avoid singularity
    bool sng = fabs(acf0[0] < 0.0001);
    if(sng) acf00 += 1.0;
    // bool sng = false;
    // if (fabs(acf0[0]) < 0.0001) {
    //   sng = true;
    //   acf0[0] += 1;
    // }

    std::copy(acf0, acf0 + ntz_, U1_);
    if(sng) U1_[0] += 1.0;
    std::fill(U1_ + ntz_, U1_ + 2 * ntz_, 0);
    vfft_->fft(U1_fft_, U1_); // U1
    std::fill(U2_, U2_ + 2 * ntz_, 0);
    std::copy(acf0 + 1, acf0 + ntz_, U2_ + 1);
    vfft_->fft(U2_fft_, U2_);  // U2

    // tr{U1'L1L1'U1}
    complex_mult(y_fft_, L1_fft_, U1_fft_, ntz_ + 1);
    vfft_->ifft(y_, y_fft_);
    trace = trace_LU(y_, y_, ntz_);

    // tr{U1'L2L2'U1}
    complex_mult(y_fft_, L2_fft_, U1_fft_, 2 * (ntz_ / 2 + 1));
    vfft_->ifft(y_, y_fft_);
    trace -= trace_LU(y_, y_, ntz_);
    // tr{U2'L1L1'U2}
    complex_mult(y_fft_, L1_fft_, U2_fft_, 2 * (ntz_ / 2 + 1));
    vfft_->ifft(y_, y_fft_);
    trace -= trace_LU(y_, y_, ntz_);
    // tr{U2'L2L2'U2}
    complex_mult(y_fft_, L2_fft_, U2_fft_, 2 * (ntz_ / 2 + 1));
    vfft_->ifft(y_, y_fft_);
    trace += trace_LU(y_, y_, ntz_);

    // trace
    trace /= gs_->delta[0];
    trace /= acf00;

    if (sng) {
      trace -= trace_inv();
      // double t0 = 0.0;
      // for (int ii = 0; ii < ntz_; ii++) {
      // 	t0 += (ntz_ - 2 * ii) * gs_->delta[ii] * gs_->delta[ii];
      // }
      // trace -= t0 / gs_->delta[0];
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
  double trace;
  double acf20 = acf2[0];
  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    if (!has_solve_) {
      solve_setup();
    }

    // bool sng = false;
    // if (fabs(acf2[0]) < 0.0001) {
    //   sng = true;
    //   acf2[0] += 1;
    // }
    bool sng = fabs(acf2[0] < 0.0001);
    if(sng) acf20 += 1.0;
		
    // we store the negative derivative of delta in vector phi_, where phi_ = solve(acf_) * toep(acf1) * delta
    product(phi_, gs_->delta, acf1);
    solve(phi_, phi_);
    trace = trace_deriv(acf2);
    if(sng) trace += trace_inv();
    trace *= -phi_[0];

    double k1;
    std::copy(phi_, phi_ + ntz_, x_);
    std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
    vfft_->fft(x_fft_, x_);
    std::copy(acf2, acf2 + ntz_, z_);
    if(sng) z_[0] += 1.0;
    std::fill(z_ + ntz_, z_ + 2 * ntz_, 0);
    vfft_->fft(z_fft_, z_);
    complex_mult(U1_fft_, x_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U1_, U1_fft_);
    std::copy(gs_->delta, gs_->delta + ntz_, y_);
    std::fill(y_ + ntz_, y_ + 2 * ntz_, 0);
    vfft_->fft(y_fft_, y_);
    complex_mult(U2_fft_, y_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U2_, U2_fft_);
    k1 = trace_LU(U1_, U2_, ntz_) / acf20;
    z_[0] = 0;
    vfft_->fft(z_fft_, z_);
    complex_mult(U1_fft_, x_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U1_, U1_fft_);
    complex_mult(U2_fft_, y_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U2_, U2_fft_);
    k1 -= trace_LU(U1_, U2_, ntz_) / acf20;

    double k2;
    x_[0] = 0;
    std::reverse(x_ + 1, x_ + ntz_);
    vfft_->fft(x_fft_, x_);
    z_[0] = acf20;
    vfft_->fft(z_fft_, z_);
    complex_mult(U1_fft_, x_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U1_, U1_fft_);
    y_[0] = 0;
    std::reverse(y_ + 1, y_ + ntz_);
    vfft_->fft(y_fft_, y_);
    complex_mult(U2_fft_, y_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U2_, U2_fft_);
    k2 = trace_LU(U1_, U2_, ntz_) / acf20;
    z_[0] = 0;
    vfft_->fft(z_fft_, z_);
    complex_mult(U1_fft_, x_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U1_, U1_fft_);
    complex_mult(U2_fft_, y_fft_, z_fft_, ntz_ + 1);
    vfft_->ifft(U2_, U2_fft_);
    k2 -= trace_LU(U1_, U2_, ntz_) / acf20;

    trace += 2 * (k1 - k2);
    trace /= gs_->delta[0];

    if (sng) {
      double t0 = 0.0;
      for (int ii = 0; ii < ntz_; ii++) {
	t0 += (ntz_ - 2 * ii) * gs_->delta[ii] * phi_[ii];
      }
      trace -= 2 * t0 / gs_->delta[0];
      trace += trace_inv() * phi_[0] / gs_->delta[0];
      // t0 = 0.0;
      // for (int ii = 0; ii < ntz_; ii++) {
      // 	t0 += (ntz_ - 2 * ii) * gs_->delta[ii] * gs_->delta[ii];
      // }
      // trace += t0 * phi_[0] / gs_->delta[0] / gs_->delta[0];
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

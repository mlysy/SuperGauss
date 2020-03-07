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
  int ntz_;        ///< Size of Toeplitz matrix.
  double* acf_;    ///< First column of the Toeplitz matrix.
  double* tzcirc_; ///< Storage for Toeplitz circulant embedding.
  double* x_;      ///< Storage for `mult/solve` input.
  double* y_;      ///< Storage for `mult/solve` output.
  bool has_acf;    ///< Boolean flag indicating whether input argument `acf_` has been modified.
  bool has_mult;   ///< Boolean flag indicating whether multiplication-related FFT has been done.
  bool has_solve;  ///< Boolean flag indicating whether inversion-related FFT has been done.
  GSchurN* Gs;     ///< Object to compute GSchur algorithm.
  VectorFFT* V1;   ///< Object for fft computations.
  
  std::complex<double>* tzcirc_fft_;
  std::complex<double>* x_fft;
  std::complex<double>* y_fft;
  double* z;
  std::complex<double>* z_fft;
  double* phi;

  std::complex<double>* L1_fft;
  std::complex<double>* L11_fft;
  std::complex<double>* L2_fft;
  std::complex<double>* L22_fft;

  // upper triangular Toeplitz matrices
  double* U1;
  double* U2;
  std::complex<double>* U1_fft;
  std::complex<double>* U2_fft;
  
  /// Prepare for multiplication.
  void mult_setup();
    
  /// Prepare for inversion.
  void solve_setup();

public:
  /// Constructor.
  Toeplitz(int n_);
  
  /// Destructor.
  ~Toeplitz();
  
  /// Setup ACF, the first column of Toeplitz matrix.
  void setAcf(const double* acfIn);
  
  /// Return ACF.
  void getAcf(double* acfOut);
  
  /// Return the size of ACF.
  int size();
  
  /// Exam whether ACF exists.
  bool hasAcf();
  
  /// Multiplication function.
  void mult(double* yOut, const double* xIn);

  /// symmetric Multiplication function.
  void product(double* yOut, const double* xIn, const double* acf1);

  /// insymmetric Multiplication function.
  void product(double* yOut, const double* xIn, const double* acf1, const double* acf2);
  
  /// Inversion function.
  void solve(double* yOut, const double* xIn);
  
  /// Determinant function.
  double logDet();

  /// Trace function.
  double Trace();

  /// Function for a trace computation
  double trace_deriv(const double* acf1);

  /// Function for another trace computation
  double trace_hess(const double* acf1, const double* acf2);

};

/// @param[in] n_ Size of Toeplitz matrix.
inline Toeplitz::Toeplitz(int n_) {
  ntz_ = n_;
  acf_ = new double[ntz_];
  has_acf = false;
  has_mult = false;
  has_solve = false;

  // GSchur algorithm only supports N > 1 case.
  if (ntz_ > 1) {
    Gs = new GSchurN(ntz_);
    V1 = new VectorFFT(2 * ntz_);

    tzcirc_ = new double[2 * ntz_]; /// storage for circulant embedding
    tzcirc_fft_ = new std::complex<double>[2 * ntz_];
    x_ = new double[2 * ntz_]; /// storage for xIn
    x_fft = new std::complex<double>[2 * ntz_];
    y_ = new double[2 * ntz_]; /// storage for yOut
    y_fft = new std::complex<double>[2 * ntz_];
    z = new double[2 * ntz_]; /// storage for temporary vector
    z_fft = new std::complex<double>[2 * ntz_];
    phi = new double[ntz_];

    L1_fft = new std::complex<double>[2 * ntz_];
    L11_fft = new std::complex<double>[2 * ntz_];
    L2_fft = new std::complex<double>[2 * ntz_];
    L22_fft = new std::complex<double>[2 * ntz_];

    U1 = new double[2 * ntz_];
    U2 = new double[2 * ntz_];
    U1_fft = new std::complex<double>[2 * ntz_];
    U2_fft = new std::complex<double>[2 * ntz_];
  }
}

inline Toeplitz::~Toeplitz() {
  delete[] acf_;

  // GSchur algorithm only supports N > 1 case.
  if (ntz_ > 1) {
    delete Gs;
    delete V1;

    delete[] tzcirc_;
    delete[] tzcirc_fft_;
    delete[] x_;
    delete[] x_fft;
    delete[] y_;
    delete[] y_fft;
    delete[] z;
    delete[] z_fft;
    delete[] phi;

    delete[] L1_fft;
    delete[] L11_fft;
    delete[] L2_fft;
    delete[] L22_fft;
    
    delete[] U1;
    delete[] U2;
    delete[] U1_fft;
    delete[] U2_fft;
  }
}


/// @param[in] acfIn Real vector of ACF.
inline void Toeplitz::setAcf(const double* acfIn) {
  std::copy(acfIn, acfIn + ntz_, acf_);
  has_acf = true;
  has_mult = false;
  has_solve = false;
  return;
}

/// @param[out] acfOut Real vector of ACF.
inline void Toeplitz::getAcf(double* acfOut) {
  std::copy(acf_, acf_ + ntz_, acfOut);
  return;
}

inline bool Toeplitz::hasAcf() { 
  return has_acf; 
}

inline int Toeplitz::size() { 
  return ntz_; 
}

/// Preparation of computing \f$T \times x\f$, i.e. generating
/// FFT([acf_, 0, rev(acf_[-1])]) and change the flag `has_mult`.
/// FFT result is stored in tzcirc_fft_.
inline void Toeplitz::mult_setup() {
  has_mult = true;
  if (ntz_ > 1) {
    std::copy(acf_, acf_ + ntz_, tzcirc_);
    tzcirc_[ntz_] = 0;
    std::copy(acf_ + 1, acf_ + ntz_, tzcirc_ + ntz_ + 1);
    std::reverse(tzcirc_ + ntz_ + 1, tzcirc_ + 2 * ntz_);
    V1->fft(tzcirc_fft_, tzcirc_);
  }
  return;
}

/// Product \f$ T \times x\f$ is computed efficiently by extending Toeplitz
/// matrix \f$T\f$ into a circulant matrix, whose multiplication with vector can
/// be done with Fast Fourier Transformation (FFT) in \f$O(n \log n) \f$ steps.
///
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
inline void Toeplitz::mult(double* yOut, const double* xIn) {
  /// tzcirc_ = [acf_, 0, rev(acf_[-1])]
  if (!has_mult) {
    mult_setup();
  }

  /// x_ = [xIn, 0]
  std::copy(xIn, xIn + ntz_, x_);
  std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
  V1->fft(x_fft, x_);

  // y_ = ifft(fft(tzcirc_) * fft(x_))[1:ntz_]
  complex_mult(y_fft, tzcirc_fft_, x_fft, 2 * (ntz_ / 2 + 1));
  V1->ifft(y_, y_fft);
  std::copy(y_, y_ + ntz_, yOut);

  return;
}



/// Preparation of computing \f$T^{-1} \times x\f$, including GSchur algorithm
/// implementation and generating FFT for Gohberg-Semencul formula.
inline void Toeplitz::solve_setup() {
  has_solve = true;

  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    Gs->compute(acf_);

    /// L11_fft stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_1'
    z[0] = Gs->delta[0];
    std::fill(z + 1, z + ntz_ + 1, 0);
    std::copy(Gs->delta + 1, Gs->delta + ntz_, z + ntz_ + 1);
    std::reverse(z + ntz_ + 1, z + 2 * ntz_);
    V1->fft(L11_fft, z);

    /// L1_fft stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_1
    std::copy(Gs->delta, Gs->delta + ntz_, z);
    std::fill(z + ntz_, z + 2 * ntz_, 0);
    V1->fft(L1_fft, z); 

    /// L22_fft stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_2'
    std::fill(z, z + ntz_ + 1, 0);
    std::copy(Gs->delta + 1, Gs->delta + ntz_, z + ntz_ + 1);
    V1->fft(L22_fft, z);
	
    /// L2_fft stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_2
    std::fill(z, z + 2 * ntz_, 0);
    std::copy(Gs->delta + 1, Gs->delta + ntz_, z + 1);
    std::reverse(z + 1, z + ntz_);
    V1->fft(L2_fft, z);

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
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
inline void Toeplitz::solve(double* yOut, const double* xIn) {
  if (!has_solve) {
    solve_setup();
  }

  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(xIn, xIn + ntz_, x_);
    std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
    V1->fft(x_fft, x_);

    // y_ = L1' * x_
    complex_mult(y_fft, L11_fft, x_fft, ntz_ + 1);
    V1->ifft(y_, y_fft);
    std::fill(y_ + ntz_, y_ + 2 * ntz_, 0);
    // z = L1 * y_ = L1 * L1' * x_
    V1->fft(y_fft, y_);
    complex_mult(z_fft, L1_fft, y_fft, ntz_ + 1);
    V1->ifft(z, z_fft);
    
    // y_ = L2' * x_
    complex_mult(y_fft, L22_fft, x_fft, ntz_ + 1);
    V1->ifft(y_, y_fft);
    std::fill(y_ + ntz_, y_ + 2 * ntz_, 0);
    // x_ = L2 * y_ = L2 * L2' * x_ (temporarily stored in x_)
    V1->fft(y_fft, y_);
    complex_mult(x_fft, L2_fft, y_fft, ntz_ + 1);
    V1->ifft(x_, x_fft);

    // yOut = (z - x_) / delta[1] = 1/delta[1] * (L1 * L1' * x_ - L2 * L2' * x_)
    for (int ii = 0; ii < ntz_; ++ii) {
      yOut[ii] = (z[ii] - x_[ii]) / Gs->delta[0];
    }
  } else {
    // N = 1 case.
    yOut[0] = xIn[0] / acf_[0];
  }
  return;
}

/// Log-determinant of Toeplitz matrix.
inline double Toeplitz::logDet() {
  if (!has_solve) {
    solve_setup();
  }
  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    return Gs->ldV;
  } else {
    // N = 1 case.
    return log(acf_[0]);
  }
}

inline double Toeplitz::Trace() {
  double t0 = 0.0;
  for (int ii = 0; ii < ntz_; ii++) {
    t0 += (ntz_ - 2 * ii) * Gs->delta[ii] * Gs->delta[ii];
  }
  return t0 / Gs->delta[0];
}

/// Computing the multiplication between a symmetric Toeplitz matrix and a vector.
/// This is a convenient function that takes outside arguments including acf1, acf2 and x and has nothing to do with the built-in argument acf_.
/// We kind of borrow the memory of Toeplitz class for such computation.
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
/// @param[in] acf1 Real vector of the first column of the asymmetric Toeplitz matrix.
inline void Toeplitz::product(double* yOut, const double* xIn, const double* acf1) {
  std::copy(acf1, acf1 + ntz_, z);
  z[ntz_] = 0;
  std::copy(acf1 + 1, acf1 + ntz_, z + ntz_ + 1);
  std::reverse(z + ntz_ + 1, z + 2 * ntz_);
  V1->fft(z_fft, z);

  /// x_ = [xIn, 0]
  std::copy(xIn, xIn + ntz_, x_);
  std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
  V1->fft(x_fft, x_);

  /// y_ = ifft(fft(z) * fft(x_))[1:ntz_]
  complex_mult(y_fft, z_fft, x_fft, ntz_ + 1);
  V1->ifft(y_, y_fft);
  std::copy(y_, y_ + ntz_, yOut);

  return;
}

/// Computing the multiplication between a non-symmetric Toeplitz matrix and a vector.
/// This is a convenient function that takes outside arguments including acf1, acf2 and x and has nothing to do with the built-in argument acf_.
/// We kind of borrow the memory of Toeplitz class for such computation.
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
/// @param[in] acf1 Real vector of the first column of the asymmetric Toeplitz matrix.
/// @param[in] acf2 Real vector of the first row of the asymmetric Toeplitz matrix.
inline void Toeplitz::product(double* yOut, const double* xIn, const double* acf1, const double* acf2) {
  /// z = [acf1, 0, rev(acf2[-1])]
  std::copy(acf1, acf1 + ntz_, z);
  z[ntz_] = 0;
  std::copy(acf2 + 1, acf2 + ntz_, z + ntz_ + 1);
  std::reverse(z + ntz_ + 1, z + 2 * ntz_);
  V1->fft(z_fft, z);

  /// x_ = [xIn, 0]
  std::copy(xIn, xIn + ntz_, x_);
  std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
  V1->fft(x_fft, x_);

  /// y_ = ifft(fft(z) * fft(x_))[1:ntz_]
  complex_mult(y_fft, z_fft, x_fft, 2 * (ntz_ / 2 + 1));
  V1->ifft(y_, y_fft);
  std::copy(y_, y_ + ntz_, yOut);

  return;
}


inline double Toeplitz::trace_deriv(const double* acf0) {
  double trace;
  double acf00 = acf0[0];

  if (ntz_ > 1) { // GSchur algorithm only supports N > 1 case.
    if (!has_solve) {
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

    std::copy(acf0, acf0 + ntz_, U1);
    if(sng) U1[0] += 1.0;
    std::fill(U1 + ntz_, U1 + 2 * ntz_, 0);
    V1->fft(U1_fft, U1); // U1
    std::fill(U2, U2 + 2 * ntz_, 0);
    std::copy(acf0 + 1, acf0 + ntz_, U2 + 1);
    V1->fft(U2_fft, U2);  // U2

    // tr{U1'L1L1'U1}
    complex_mult(y_fft, L1_fft, U1_fft, ntz_ + 1);
    V1->ifft(y_, y_fft);
    trace = trace_LU(y_, y_, ntz_);

    // tr{U1'L2L2'U1}
    complex_mult(y_fft, L2_fft, U1_fft, 2 * (ntz_ / 2 + 1));
    V1->ifft(y_, y_fft);
    trace -= trace_LU(y_, y_, ntz_);
    // tr{U2'L1L1'U2}
    complex_mult(y_fft, L1_fft, U2_fft, 2 * (ntz_ / 2 + 1));
    V1->ifft(y_, y_fft);
    trace -= trace_LU(y_, y_, ntz_);
    // tr{U2'L2L2'U2}
    complex_mult(y_fft, L2_fft, U2_fft, 2 * (ntz_ / 2 + 1));
    V1->ifft(y_, y_fft);
    trace += trace_LU(y_, y_, ntz_);

    // trace
    trace /= Gs->delta[0];
    trace /= acf00;

    if (sng) {
      trace -= Trace();
      // double t0 = 0.0;
      // for (int ii = 0; ii < ntz_; ii++) {
      // 	t0 += (ntz_ - 2 * ii) * Gs->delta[ii] * Gs->delta[ii];
      // }
      // trace -= t0 / Gs->delta[0];
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
  double trace2;
  double acf20 = acf2[0];
  if (ntz_ > 1) {
    // GSchur algorithm only supports N > 1 case.
    if (!has_solve) {
      solve_setup();
    }

    // bool sng = false;
    // if (fabs(acf2[0]) < 0.0001) {
    //   sng = true;
    //   acf2[0] += 1;
    // }
    bool sng = fabs(acf2[0] < 0.0001);
    if(sng) acf20 += 1.0;
		
    // we store the negative derivative of delta in vector phi, where phi = solve(acf_) * toep(acf1) * delta
    product(phi, Gs->delta, acf1);
    solve(phi, phi);
    trace2 = trace_deriv(acf2);
    if(sng) trace2 += Trace();
    trace2 *= -phi[0];

    double k1;
    std::copy(phi, phi + ntz_, x_);
    std::fill(x_ + ntz_, x_ + 2 * ntz_, 0);
    V1->fft(x_fft, x_);
    std::copy(acf2, acf2 + ntz_, z);
    if(sng) z[0] += 1.0;
    std::fill(z + ntz_, z + 2 * ntz_, 0);
    V1->fft(z_fft, z);
    complex_mult(U1_fft, x_fft, z_fft, ntz_ + 1);
    V1->ifft(U1, U1_fft);
    std::copy(Gs->delta, Gs->delta + ntz_, y_);
    std::fill(y_ + ntz_, y_ + 2 * ntz_, 0);
    V1->fft(y_fft, y_);
    complex_mult(U2_fft, y_fft, z_fft, ntz_ + 1);
    V1->ifft(U2, U2_fft);
    k1 = trace_LU(U1, U2, ntz_) / acf20;
    z[0] = 0;
    V1->fft(z_fft, z);
    complex_mult(U1_fft, x_fft, z_fft, ntz_ + 1);
    V1->ifft(U1, U1_fft);
    complex_mult(U2_fft, y_fft, z_fft, ntz_ + 1);
    V1->ifft(U2, U2_fft);
    k1 -= trace_LU(U1, U2, ntz_) / acf20;

    double k2;
    x_[0] = 0;
    std::reverse(x_ + 1, x_ + ntz_);
    V1->fft(x_fft, x_);
    z[0] = acf20;
    V1->fft(z_fft, z);
    complex_mult(U1_fft, x_fft, z_fft, ntz_ + 1);
    V1->ifft(U1, U1_fft);
    y_[0] = 0;
    std::reverse(y_ + 1, y_ + ntz_);
    V1->fft(y_fft, y_);
    complex_mult(U2_fft, y_fft, z_fft, ntz_ + 1);
    V1->ifft(U2, U2_fft);
    k2 = trace_LU(U1, U2, ntz_) / acf20;
    z[0] = 0;
    V1->fft(z_fft, z);
    complex_mult(U1_fft, x_fft, z_fft, ntz_ + 1);
    V1->ifft(U1, U1_fft);
    complex_mult(U2_fft, y_fft, z_fft, ntz_ + 1);
    V1->ifft(U2, U2_fft);
    k2 -= trace_LU(U1, U2, ntz_) / acf20;

    trace2 += 2 * (k1 - k2);
    trace2 /= Gs->delta[0];

    if (sng) {
      double t0 = 0.0;
      for (int ii = 0; ii < ntz_; ii++) {
	t0 += (ntz_ - 2 * ii) * Gs->delta[ii] * phi[ii];
      }
      trace2 -= 2 * t0 / Gs->delta[0];
      t0 = 0.0;
      for (int ii = 0; ii < ntz_; ii++) {
	t0 += (ntz_ - 2 * ii) * Gs->delta[ii] * Gs->delta[ii];
      }
      trace2 += t0 * phi[0] / Gs->delta[0] / Gs->delta[0];
      // acf2[0] -= 1;
    }
		
  }
  else {
    // N = 1 case.
    trace2 = acf1[0] * acf2[0] / acf_[0] / acf_[0];
  }
  return trace2;
}


#endif

/// @file Toeplitz.h

///////////////////////////////////////////////
// Toeplitz matrix class
///////////////////////////////////////////////

#ifndef Toeplitz_h
#define Toeplitz_h 1

#include "GSchur.h"

/// Class for Toeplitz-related computation.
///
/// Generate the memory for Toeplitz computation including multiplication,
/// inversion, determinant, first and second order derivative of the
/// log-determinant.
class Toeplitz {
  int n;               ///< Integer indicating the size of Toeplitz matrix.
  VectorFFT* L1fft;    ///< VectorFFT space for matrix \f$L_1\f$ for inversion.
  VectorFFT* L11fft;   ///< VectorFFT space for matrix \f$L_1'\f$ for inversion.
  VectorFFT* L2fft;    ///< VectorFFT space for matrix \f$L_2\f$ for inversion.
  VectorFFT* L22fft;   ///< VectorFFT space for matrix \f$L_2'\f$ for inversion.
  VectorFFT* xfft;     ///< VectorFFT space for vector x.
  VectorFFT* Lxfft;    ///< VectorFFT space for vector \f$L_1' x\f$.
  VectorFFT* U1fft;    ///< VectorFFT space for matrix \f$U_1\f$.
  VectorFFT* U2fft;    ///< VectorFFT space for matrix \f$U_2\f$.
  VectorFFT* Toepfft;  ///< VectorFFT space for matrix \f$T\f$.
  VectorIFFT* Invfft;  ///< VectorIFFT space for multiple usage.
  double* phi2;        ///< Vector for second derivative of log-determinant.
  double* temVec;      ///< Vector for second derivative of log-determinant.
  bool has_acf;    ///< Boolen, flag indicating whether input argument `acf` has
                   ///< been modified.
  bool has_mult;   ///< Boolen, flag indicating whether multiplication-related
                   ///< FFT has been done.
  bool has_solve;  ///< Boolen, flag indicating whether inversion-related FFT
                   ///< has been done.
  GSchurN* Gs;     ///< GSchur Class.
  double* acf;     ///< Vector of input, the first column of Toeplitz matrix.

  /// Check whether vector is all-zero.
  bool all_zeros(double*, int);
  /// Compute the trace of product of lower and upper triangular matrices.
  double trace_LU(double*, double*, int);
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
  void setAcf(double* acfIn);
  /// Return ACF.
  void getAcf(double* acfOut);
  /// Return the size of ACF.
  int size();
  /// Exam whether ACF exists.
  bool hasAcf();
  /// Multiplication function.
  void multVec(double* yOut, double* xIn);
  /// Inversion function.
  void solveVec(double* yOut, double* xIn);
  /// Determinant function.
  double logDet();
  /// Trace function for derivative of log-determinant.
  double traceProd(double*);
  /// Trace function for second derivative of log-determinant.
  double traceDeriv(double*, double*);
  /// Obtain the result of GSchur algorithm.
  void getPhi(double* phiOut);
};

/// @param[in] n_ Size of Toeplitz matrix.
inline Toeplitz::Toeplitz(int n_) {
  n = n_;
  acf = new double[n];
  phi2 = new double[n];
  temVec = new double[n];
  has_acf = false;

  // GSchur algorithm only supports N > 1 case.
  if (n_ > 1) {
    Gs = new GSchurN(n, 64);  // default binary modulus 64
    L1fft = new VectorFFT(2 * n);
    L11fft = new VectorFFT(2 * n);
    L2fft = new VectorFFT(2 * n);
    L22fft = new VectorFFT(2 * n);
    xfft = new VectorFFT(2 * n);
    Lxfft = new VectorFFT(2 * n);
    U1fft = new VectorFFT(2 * n);
    U2fft = new VectorFFT(2 * n);
    Toepfft = new VectorFFT(2 * n);
    Invfft = new VectorIFFT(2 * n);
  }
}

inline Toeplitz::~Toeplitz() {
  delete[] acf;
  delete[] phi2;
  delete[] temVec;

  // GSchur algorithm only supports N > 1 case.
  if (n > 1) {
    delete Gs;
    delete L1fft;
    delete L11fft;
    delete L2fft;
    delete L22fft;
    delete xfft;
    delete Lxfft;
    delete U1fft;
    delete U2fft;
    delete Toepfft;
    delete Invfft;
  }
}

/// Compute \f$\text{tr}(L U)\f$, where \f$L\f$ is a lower triangular Toeplitz
/// matrix and \f$U\f$ is a upper triangular Toeplitz matrix. It is an
/// \f$O(n)\f$ algorithm.
///
/// @param[in] acf1 Real vector of length `n`, first column of \f$L\f$.
/// @param[in] acf2 Real vector of length `n`, first row of \f$U\f$.
/// @param[in] n Integer indicating the length of input.
inline double Toeplitz::trace_LU(double* acf1, double* acf2, int n) {
  double trace = 0;
  for (int ii = 0; ii < n; ++ii) {
    trace += (n - ii) * acf1[ii] * acf2[ii];
  }
  return trace;
}

/// @param[in] acf Real input vector.
/// @param[in] n Integer for length of input vector.
/// @return Boolen indicating whether vector is all zero or not.
inline bool Toeplitz::all_zeros(double* acf, int n) {
  bool acf_is_0 = true;
  for (int ii = 0; ii < n; ++ii) {
    if (fabs(acf[ii]) > 0.0) {
      acf_is_0 = false;
      break;
    }
  }
  return acf_is_0;
}

/// @param[in] acfIn Real vector of ACF.
inline void Toeplitz::setAcf(double* acfIn) {
  std::copy(acfIn, acfIn + n, acf);
  has_acf = true;
  has_mult = false;
  has_solve = false;
  return;
}

/// @param[out] acfOut Real vector of ACF.
inline void Toeplitz::getAcf(double* acfOut) {
  std::copy(acf, acf + n, acfOut);
  return;
}

inline bool Toeplitz::hasAcf() { return has_acf; }

inline int Toeplitz::size() { return n; }

/// GSchur algorithm returns the first column of inverse Toeplitz matrix:
/// \f$\rho\f$.
///
/// @param[out] phiOut Real vector
inline void Toeplitz::getPhi(double* phiOut) {
  std::copy(Gs->Phi, Gs->Phi + n, phiOut);
  return;
}

/// Preparation of computing \f$T \times x\f$, including pasting, generating
/// FFT(acf) and change the flag `has_mult`.
inline void Toeplitz::mult_setup() {
  has_mult = true;
  if (n > 1) {
    std::copy(acf, acf + n, Toepfft->in);
    std::copy(acf + 1, acf + n, Toepfft->in + n + 1);
    std::reverse(Toepfft->in + n + 1, Toepfft->in + 2 * n);
    Toepfft->fft();
  }
  return;
}

/// Product \f$ T \times x\f$ is computed efficiently by extending Toeplitz
/// matrix \f$T\f$ into a circulant matrix, whose multiplication with vector can
/// be done with Fast Fourier Transformation (FFT) in \f$O(n \log n) \f$ steps.
///
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
inline void Toeplitz::multVec(double* yOut, double* xIn) {
  bool acf_is_0 = all_zeros(acf, n);
  if (acf_is_0) {
    // zero matrix times vector
    std::fill(yOut, yOut + n, 0);
    return;
  }
  if (!has_mult) {
    mult_setup();
  }
  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(xIn, xIn + n, xfft->in);
    xfft->fft();
    vecConv(Invfft->in, Toepfft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, yOut);
  } else {
    // N = 1 case.
    yOut[0] = acf[0] * xIn[0];
  }
  return;
}

/// Preparation of computing \f$T^{-1} \times x\f$, including GSchur algorithm
/// implementation and generating FFT for Gohberg-Semencul formula.
inline void Toeplitz::solve_setup() {
  has_solve = true;
  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    Gs->Compute(acf);
    L11fft->in[0] = Gs->Phi[0];
    std::copy(Gs->Phi + 1, Gs->Phi + n, L11fft->in + n + 1);
    std::reverse(L11fft->in + n + 1, L11fft->in + 2 * n);
    L11fft->fft();  // L1'
    std::copy(Gs->Phi, Gs->Phi + n, L1fft->in);
    L1fft->fft();  // L1
    std::copy(Gs->Phi + 1, Gs->Phi + n, L22fft->in + n + 1);
    L22fft->fft();  // L2'
    std::copy(Gs->Phi + 1, Gs->Phi + n, L2fft->in + 1);
    std::reverse(L2fft->in + 1, L2fft->in + n);
    L2fft->fft();  // L2
  }
  return;
}

/// Multiplication between inverse Toeplitz matrix and vector: \f$T^{-1} \times
/// x\f$ involves the Gohberg-Semencul formula which decomposes inverse Toeplitz
/// into sum of products of lower and upper triangular Toeplitz matrices: \f$
/// T^{-1} = \frac{1}{\rho_0} [L_1'L_1 - L_2' L_2] \f$, where \f$L_1, L_2\f$ are
/// lower-trangular Toeplitz matrices and \f$\rho_0\f$ is the first element of
/// GSchur Algorithm output: length-N vector \f$\rho\f$.
///
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
inline void Toeplitz::solveVec(double* yOut, double* xIn) {
  // if (!has_solve) {
    solve_setup();
  // }

  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(xIn, xIn + n, xfft->in);
    xfft->fft();

    // L1 * L1' * x
    vecConv(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
    Lxfft->fft();
    vecConv(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, yOut);  // stored in yOut

    // L2 * L2' * x
    vecConv(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
    Lxfft->fft();
    vecConv(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    // 1/sigma2 (L1 L1'x - L2 L2'x)
    for (int ii = 0; ii < n; ++ii) {
      yOut[ii] -= Invfft->out[ii];
      yOut[ii] /= Gs->Phi[0];
    }
  } else {
    // N = 1 case.
    yOut[0] = xIn[0] / acf[0];
  }
  return;
}

/// Log-determinant of Toeplitz matrix.
inline double Toeplitz::logDet() {
  if (!has_solve) {
    solve_setup();
  }
  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    return Gs->ldV;
  } else {
    // N = 1 case.
    return log(acf[0]);
  }
}

/// First order derivative of log-determinant is:
/// \f[
///   \frac{\partial \log\mid T^{-1}\mid}{\partial \theta} = \text{tr} (T^{-1}
///   V), \text{  where  } V = \frac{\partial T}{\partial \theta}
/// \f]
/// V is also Toeplitz since it is the derivative of Toeplitz matrix, and has
/// displacement rank 2 :\f$V = \frac{1}{t_0} [U_1 U_1' - U_2 U_2']\f$, where
/// \f$U_1\f$ is a upper triangular Toeplitz matrix. Thus \f{eqnarray*}{
///   \text{tr} (T^{-1} V) &=& \frac{1}{\rho_0 t_0} \text{tr} \{ [L_1'L_1 - L_2'
///   L_2][U_1 U_1' - U_2 U_2'] \} \\
///                        &=& \frac{1}{\rho_0 t_0} (\text{tr} \{ U_1'L_1 L_1'
///                        U_1 \} + ...) \\
///                        &=& \frac{1}{\rho_0 t_0} (\text{tr} \{ A_1 A_1' \} +
///                        ...)
/// \f}
/// where \f$A = U_1'L_1\f$ is the product of two lower triangular Toeplitz
/// matrices and is also a lower triangular Toeplitz matrix, making computing
/// \f$\text{tr} \{ A_1 A_1' \}\f$ efficient.
///
/// If the first element of Toeplitz matrix \f$V\f$: \f$V[1,1]\f$ is very
/// small(criteria is `0.0001` in coding), we use following steps to stabilize
/// the algorithm:
/// \f[
///   \text{tr}(T^{-1} V) = \text{tr}(T^{-1} (V+I)) - \text{tr}(T^{-1})
/// \f]
/// where \f$I\f$ is an identity matrix.
///
/// @param[in] acf0 Real vector of first column of Toeplitz matrix \f$V\f$.
/// @return Real number of trace.
inline double Toeplitz::traceProd(double* acf0) {
  double trace;
  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    bool is_small_acf0;
    double t0;
    if (all_zeros(acf0, n)) {
      return 0.0;
    }
    if (!has_solve) {
      solve_setup();
    }
    std::copy(acf0 + 1, acf0 + n, U1fft->in + 1);
    is_small_acf0 = fabs(acf0[0]) < 0.0001;
    t0 = is_small_acf0 ? acf0[0] + 1.0 : acf0[0];
    U1fft->in[0] = t0;
    U1fft->fft();  // U1
    std::copy(acf0 + 1, acf0 + n, U2fft->in + 1);
    U2fft->fft();  // U2

    // tr{U1'L1L1'U1}
    vecConv(Invfft->in, L1fft->out, U1fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace = trace_LU(Invfft->out, Invfft->out, n);
    // tr{U1'L2L2'U1}
    vecConv(Invfft->in, L2fft->out, U1fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace -= trace_LU(Invfft->out, Invfft->out, n);
    // tr{U2'L1L1'U2}
    vecConv(Invfft->in, L1fft->out, U2fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace -= trace_LU(Invfft->out, Invfft->out, n);
    // tr{U2'L2L2'U2}
    vecConv(Invfft->in, L2fft->out, U2fft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace += trace_LU(Invfft->out, Invfft->out, n);

    // trace
    trace /= Gs->Phi[0];
    trace /= t0;
    if (is_small_acf0) {
      t0 = 0.0;
      for (int ii = 1; ii < n; ii++) {
        t0 += (n - 2 * ii) * Gs->Phi[ii] * Gs->Phi[ii];
      }
      trace -= n * Gs->Phi[0] + t0 / Gs->Phi[0];
    }
  } else {
    // N = 1 case.
    trace = acf0[0] / acf[0];
  }
  return trace;
}

/// Second order derivative of log-determinant is:
/// \f[
///   \frac{\partial^2 \log\mid T^{-1}\mid}{\partial \theta_1 \partial \theta_2}
///   = \text{tr} (T^{-1} V_0) - \text{tr} (T^{-1} V_1 T^{-1} V_2), \text{ where
///   } V_0 = \frac{\partial^2 T}{\partial \theta_1 \partial \theta_2}, V_1 =
///   \frac{\partial T}{\partial \theta_1}, V_2 = \frac{\partial T}{\partial
///   \theta_2}
/// \f]
/// \f$V_0, V_1, V_2\f$ are Toeplitz matrices. Computation of \f$\text{tr}
/// (T^{-1} V_1 T^{-1} V_2)\f$ is non-trival and details are in paper.
///
/// If the first element of Toeplitz matrix \f$V_2\f$: \f$V_2[1,1]\f$ is very
/// small(criteria is `0.0001` in coding), we use following steps
/// to stabilize the algorithm:
///
/// \f[
///   \text{tr}(T^{-1} V_1 T^{-1} V_2) = \text{tr}(T^{-1} V_1 T^{-1} (V_2+I)) -
///   \text{tr}(T^{-1} V_1 T^{-1})
/// \f]
/// where \f$I\f$ is an identity matrix.
///
/// @param[in] acf1 Real vector of first column of Toeplitz matrix \f$V_1\f$.
/// @param[in] acf2 Real vector of first column of Toeplitz matrix \f$V_2\f$.
/// @return Real number of trace.
inline double Toeplitz::traceDeriv(double* acf1, double* acf2) {
  double trace2;
  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    double t0;
    bool is_small_acf2;
    if (all_zeros(acf1, n) | all_zeros(acf2, n)) {
      return 0.0;
    }
    if (!has_solve) {
      solve_setup();
    }
    // check first term
    is_small_acf2 = fabs(acf2[0]) < 0.0001;
    t0 = is_small_acf2 ? acf2[0] + 1.0 : acf2[0];

    // (1) phi2 = - T^-1 * V_1 * phi:

    // phi2 = V_1 * phi
    std::copy(acf1, acf1 + n, Lxfft->in);
    std::copy(acf1 + 1, acf1 + n, Lxfft->in + n + 1);
    std::reverse(Lxfft->in + n + 1, Lxfft->in + 2 * n);
    Lxfft->fft();
    std::fill(Lxfft->in + n, Lxfft->in + 2 * n,
              0);  // (readjust back to 0 assignment for latter half)
    std::copy(Gs->Phi, Gs->Phi + n, xfft->in);
    xfft->fft();
    vecConv(Invfft->in, Lxfft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, phi2);
    // phi2 = - T^-1 * phi2
    std::copy(phi2, phi2 + n, xfft->in);
    xfft->fft();

    // L1 * L1' * x
    vecConv(Invfft->in, L11fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
    Lxfft->fft();
    vecConv(Invfft->in, L1fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, phi2);

    // L2 * L2' * x
    vecConv(Invfft->in, L22fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, Lxfft->in);
    Lxfft->fft();
    vecConv(Invfft->in, L2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();

    // phi2 = -1/phi[1] * (L1 L1'x - L2 L2'x)
    for (int ii = 0; ii < n; ++ii) {
      phi2[ii] -= Invfft->out[ii];
      phi2[ii] /= Gs->Phi[0];
    }

    // tr = phi2[1] * tr(T^-1 * V_2) + 2 * tr(L1(-phi2) * L1(phi)' * V_2) - 2 *
    // tr(_L2(-phi2) * _L2(phi)' * V_2)

    // (1) phi2[1] * tr(T^-1 * V_2)
    trace2 = phi2[0] * traceProd(acf2);

    // (2) tr(L1(-phi2) * L1(phi)' * V_2) = -tr(L1(phi2) * L1(phi)' * L1'(acf2)
    // * L1(acf2)) + tr(L1(phi2) * L1(phi)' * L2'(acf2) * L2(acf2))

    // (2.1) tr(L1(phi2) * L1(phi)' * L1'(acf2) * L1(acf2)) = trace_LU(L1(acf2)
    // * L1(phi2), L1(acf2) * L1(phi)) = temVec
    std::copy(acf2, acf2 + n, U1fft->in);
    U1fft->in[0] = t0;
    U1fft->fft();
    std::copy(phi2, phi2 + n, Lxfft->in);
    Lxfft->fft();
    vecConv(Invfft->in, U1fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, temVec);
    std::copy(Gs->Phi, Gs->Phi + n, xfft->in);
    xfft->fft();
    vecConv(Invfft->in, U1fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace2 -= 2 * trace_LU(Invfft->out, temVec, n) / t0;

    // (2.2) tr(L1(phi2) * L1(phi)' * L2'(acf2) * L2(acf2)) = trace_LU(L2(acf2)
    // * L1(phi2), L2(acf2) * L1(phi)) = temVec
    std::copy(acf2 + 1, acf2 + n, U2fft->in + 1);
    U2fft->in[0] = 0;
    U2fft->fft();
    vecConv(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, temVec);
    vecConv(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace2 += 2 * trace_LU(Invfft->out, temVec, n) / t0;

    // (3) tr(_L2(phi2) * _L2(phi)' * Toeplitz_3) = tr(_L2(phi2) * _L2'(phi) *
    // L1'(acf2) * L1(acf2)) / acf2[1] - tr(_L2(phi2) * _L2(phi)' * L2'(acf2) *
    // L2(acf2)) / acf2[1]

    // (3.1) tr(L2(phi2) * L2'(phi) * L1'(acf2) * L1(acf2)) = trace_LU(L1(acf2)
    // * _L2(phi2), L1(acf2) * _L2(phi)) = temVec
    Lxfft->in[0] = 0;
    std::copy(phi2 + 1, phi2 + n, Lxfft->in + 1);
    std::reverse(Lxfft->in + 1, Lxfft->in + n);
    Lxfft->fft();
    vecConv(Invfft->in, U1fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, temVec);
    xfft->in[0] = 0;
    std::copy(Gs->Phi + 1, Gs->Phi + n, xfft->in + 1);
    std::reverse(xfft->in + 1, xfft->in + n);
    xfft->fft();
    vecConv(Invfft->in, U1fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace2 += 2 * trace_LU(Invfft->out, temVec, n) / t0;

    // (3.2) tr(_L2(phi2) * _L2(phi)' * L2'(acf2) * L2(acf2)) =
    // trace_LU(L2(acf2) * _L2(phi2), L2(acf2) * _L2(phi)) = temVec
    vecConv(Invfft->in, U2fft->out, Lxfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    std::copy(Invfft->out, Invfft->out + n, temVec);
    vecConv(Invfft->in, U2fft->out, xfft->out, 2 * (n / 2 + 1));
    Invfft->Ifft();
    trace2 -= 2 * trace_LU(Invfft->out, temVec, n) / t0;

    // trace2 = trace2 / phi[1]
    trace2 /= Gs->Phi[0];
    trace2 *= -1.0;
    if (is_small_acf2) {
      t0 = 0.0;
      for (int ii = 1; ii < n; ii++) {
        t0 += (n - 2 * ii) * Gs->Phi[ii] * phi2[ii];
      }
      trace2 -= 2 * n * phi2[0] + 2 * t0 / Gs->Phi[0];
    }
  } else {
    // N = 1 case.
    trace2 = acf1[0] * acf2[0] / acf[0] / acf[0];
  }
  return trace2;
}

#endif

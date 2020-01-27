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
protected:
  int n;               ///< Integer indicating the size of Toeplitz matrix.
  
  bool has_acf;    ///< Boolen, flag indicating whether input argument `acf` has been modified.
  bool has_mult;   ///< Boolen, flag indicating whether multiplication-related FFT has been done.
  bool has_solve;  ///< Boolen, flag indicating whether inversion-related FFT has been done.

  GSchurN* Gs;     ///< GSchur Class.
  VectorFFT* V1;   ///< VectorFFT Class for fft computations.
  double* acf;     ///< Vector of input, the first column of Toeplitz matrix.
  
  double* Tz;
  std::complex<double>* Tz_fft;
  double* x;
  std::complex<double>* x_fft;
  double* y;
  std::complex<double>* y_fft;
  double* z;
  std::complex<double>* z_fft;

  std::complex<double>* L1_fft;
  std::complex<double>* L11_fft;
  std::complex<double>* L2_fft;
  std::complex<double>* L22_fft;
  
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
  void mult(double* yOut, double* xIn);

  /// symmetric Multiplication function.
  void product(double* yOut, double* xIn, double* acf1);

  /// insymmetric Multiplication function.
  void product(double* yOut, double* xIn, double* acf1, double* acf2);
  
  /// Inversion function.
  void solve(double* yOut, double* xIn);
  
  /// Determinant function.
  double logDet();
};

/// @param[in] n_ Size of Toeplitz matrix.
inline Toeplitz::Toeplitz(int n_) {
  n = n_;
  acf = new double[n];
  has_acf = false;
  has_mult = false;
  has_solve = false;

  // GSchur algorithm only supports N > 1 case.
  if (n_ > 1) {
    Gs = new GSchurN(n);
	V1 = new VectorFFT(2 * n);

	Tz = new double[2 * n]; /// storage for circulant embedding
	Tz_fft = new std::complex<double>[2 * n];
	x = new double[2 * n]; /// storage for xIn
	x_fft = new std::complex<double>[2 * n];
	y = new double[2 * n]; /// storage for yOut
	y_fft = new std::complex<double>[2 * n];
	z = new double[2 * n]; /// storage for temporary vector
	z_fft = new std::complex<double>[2 * n];

	L1_fft = new std::complex<double>[2 * n];
	L11_fft = new std::complex<double>[2 * n];
	L2_fft = new std::complex<double>[2 * n];
	L22_fft = new std::complex<double>[2 * n];

  }
}

inline Toeplitz::~Toeplitz() {
  delete[] acf;

  // GSchur algorithm only supports N > 1 case.
  if (n > 1) {
    delete Gs;
	delete V1;

	delete[] Tz;
	delete[] Tz_fft;
	delete[] x;
	delete[] x_fft;
	delete[] y;
	delete[] y_fft;
	delete[] z;
	delete[] z_fft;

	delete[] L1_fft;
	delete[] L11_fft;
	delete[] L2_fft;
	delete[] L22_fft;
  }
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

inline bool Toeplitz::hasAcf() { 
	return has_acf; 
}

inline int Toeplitz::size() { 
	return n; 
}

/// Preparation of computing \f$T \times x\f$, i.e. generating
/// FFT([acf, 0, rev(acf[-1])]) and change the flag `has_mult`.
/// FFT result is stored in Tz_fft.
inline void Toeplitz::mult_setup() {
  has_mult = true;
  if (n > 1) {
	std::copy(acf, acf + n, Tz);
	Tz[n] = 0;
    std::copy(acf + 1, acf + n, Tz + n + 1);
    std::reverse(Tz + n + 1, Tz + 2 * n);
	V1->fft(Tz_fft, Tz);
  }
  return;
}

/// Product \f$ T \times x\f$ is computed efficiently by extending Toeplitz
/// matrix \f$T\f$ into a circulant matrix, whose multiplication with vector can
/// be done with Fast Fourier Transformation (FFT) in \f$O(n \log n) \f$ steps.
///
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
inline void Toeplitz::mult(double* yOut, double* xIn) {
	/// Tz = [acf, 0, rev(acf[-1])]
	if (!has_mult) {
		mult_setup();
	}

	/// x = [xIn, 0]
	std::copy(xIn, xIn + n, x);
	std::fill(x + n, x + 2 * n, 0);
	V1->fft(x_fft, x);

	// y = ifft(fft(Tz) * fft(x))[1:n]
	Complex_Mult(y_fft, Tz_fft, x_fft, 2 * (n / 2 + 1));
	V1->ifft(y, y_fft);
	std::copy(y, y + n, yOut);

	return;
}



/// Preparation of computing \f$T^{-1} \times x\f$, including GSchur algorithm
/// implementation and generating FFT for Gohberg-Semencul formula.
inline void Toeplitz::solve_setup() {
  has_solve = true;

  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    Gs->Compute(acf);

	/// L11_fft stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_1'
	z[0] = Gs->delta[0];
	std::fill(z + 1, z + n + 1, 0);
	std::copy(Gs->delta + 1, Gs->delta + n, z + n + 1);
	std::reverse(z + n + 1, z + 2 * n);
	V1->fft(L11_fft, z);

	/// L1_fft stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_1
	std::copy(Gs->delta, Gs->delta + n, z);
	std::fill(z + n, z + 2 * n, 0);
	V1->fft(L1_fft, z); 

	/// L22_fft stores the fft of the first column of the circulant embedding of upper triangular Toeplitz matrix L_2'
	std::fill(z, z + n + 1, 0);
	std::copy(Gs->delta + 1, Gs->delta + n, z + n + 1);
	V1->fft(L22_fft, z);
	
	/// L2_fft stores the fft of the first column of the circulant embedding of lower triangular Toeplitz matrix L_2
	std::fill(z, z + 2 * n, 0);
	std::copy(Gs->delta + 1, Gs->delta + n, z + 1);
	std::reverse(z + 1, z + n);
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
inline void Toeplitz::solve(double* yOut, double* xIn) {
  if (!has_solve) {
    solve_setup();
  }

  if (n > 1) {
    // GSchur algorithm only supports N > 1 case.
    std::copy(xIn, xIn + n, x);
	std::fill(x + n, x + 2 * n, 0);
	V1->fft(x_fft, x);

    // y = L1' * x
	Complex_Mult(y_fft, L11_fft, x_fft, n + 1);
	V1->ifft(y, y_fft);
	std::fill(y + n, y + 2 * n, 0);
	// z = L1 * y = L1 * L1' * x
	V1->fft(y_fft, y);
	Complex_Mult(z_fft, L1_fft, y_fft, n + 1);
	V1->ifft(z, z_fft);
    
    // y = L2' * x
	Complex_Mult(y_fft, L22_fft, x_fft, n + 1);
	V1->ifft(y, y_fft);
	std::fill(y + n, y + 2 * n, 0);
	// x = L2 * y = L2 * L2' * x (temporarily stored in x)
	V1->fft(y_fft, y);
	Complex_Mult(x_fft, L2_fft, y_fft, n + 1);
	V1->ifft(x, x_fft);

    // yOut = (z - x) / delta[1] = 1/delta[1] * (L1 * L1' * x - L2 * L2' * x)
    for (int ii = 0; ii < n; ++ii) {
      yOut[ii] = (z[ii] - x[ii]) / Gs->delta[0];
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

/// Computing the multiplication between a symmetric Toeplitz matrix and a vector.
/// This is a convenient function that takes outside arguments including acf1, acf2 and x and has nothing to do with the built-in argument acf.
/// We kind of borrow the memory of Toeplitz class for such computation.
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
/// @param[in] acf1 Real vector of the first column of the asymmetric Toeplitz matrix.
inline void Toeplitz::product(double* yOut, double* xIn, double* acf1) {
	std::copy(acf1, acf1 + n, z);
	z[n] = 0;
	std::copy(acf1 + 1, acf1 + n, z + n + 1);
	std::reverse(z + n + 1, z + 2 * n);
	V1->fft(z_fft, z);

	/// x = [xIn, 0]
	std::copy(xIn, xIn + n, x);
	std::fill(x + n, x + 2 * n, 0);
	V1->fft(x_fft, x);

	/// y = ifft(fft(z) * fft(x))[1:n]
	Complex_Mult(y_fft, z_fft, x_fft, n + 1);
	V1->ifft(y, y_fft);
	std::copy(y, y + n, yOut);

	return;
}

/// Computing the multiplication between a non-symmetric Toeplitz matrix and a vector.
/// This is a convenient function that takes outside arguments including acf1, acf2 and x and has nothing to do with the built-in argument acf.
/// We kind of borrow the memory of Toeplitz class for such computation.
/// @param[out] yOut Real vector of result.
/// @param[in] xIn Real vector for product.
/// @param[in] acf1 Real vector of the first column of the asymmetric Toeplitz matrix.
/// @param[in] acf2 Real vector of the first row of the asymmetric Toeplitz matrix.
inline void Toeplitz::product(double* yOut, double* xIn, double* acf1, double* acf2) {
	/// z = [acf1, 0, rev(acf2[-1])]
	std::copy(acf1, acf1 + n, z);
	z[n] = 0;
	std::copy(acf2 + 1, acf2 + n, z + n + 1);
	std::reverse(z + n + 1, z + 2 * n);
	V1->fft(z_fft, z);

	/// x = [xIn, 0]
	std::copy(xIn, xIn + n, x);
	std::fill(x + n, x + 2 * n, 0);
	V1->fft(x_fft, x);

	/// y = ifft(fft(z) * fft(x))[1:n]
	Complex_Mult(y_fft, z_fft, x_fft, 2 * (n / 2 + 1));
	V1->ifft(y, y_fft);
	std::copy(y, y + n, yOut);

	return;
}

#endif

/// @file VectorFFT.h

///////////////////////////////////////////////
// Fast Fourier Transformation and Inverse
///////////////////////////////////////////////

#ifndef VF_h
#define VF_h 1

// usual header
#include <Rcpp.h>
#include <complex>
#include <fftw3.h>
#include <iostream>

/// Forward FFT from real to complex.
class VF {
 private:
  fftw_plan planfor;  ///< plan for forward FFT.
  fftw_plan planback;  ///< plan for backward FFT.

  double* x_;
  fftw_complex* y_;
  int nr;
 public:
  int n_;         ///< Size of input vector.
  
  /// Constructor.
  VF(int);
  /// Destructor.
  ~VF();
  /// Perform the FFT on the input data.
  void fft(std::complex<double>*, double*);
  /// Perform the inverse FFT on the input data.
  void ifft(double*, std::complex<double>*);
  /// Return the pointer of x_
  double* get_ptr_x_();
  /// Return the pointer of y_
  fftw_complex* get_ptr_y_();

};

/// @param[in] n Size of the input vector.
inline VF::VF(int n) {
  n_ = n;
  nr = ceil((double)(n + 1) / 2);
  x_ = new double[n];
  std::fill(x_, x_ + n, 0);
  y_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  planfor = fftw_plan_dft_r2c_1d(n, x_, y_, FFTW_ESTIMATE);
  planback = fftw_plan_dft_c2r_1d(n, y_, x_, FFTW_ESTIMATE);
  return;
}

inline VF::~VF() {
  delete[] x_;
  fftw_free(y_);
  fftw_destroy_plan(planfor);
  fftw_destroy_plan(planback);
}

/// The result gets stored in the object member `out` of type `fftw_complex`,
/// which is a 2-d array.  Elements are accessed via e.g., `out[0][1]`,
/// `out[n-1][0]`.
inline void VF::fft(std::complex<double>* y, double* x) {
	std::copy(x, x + n_, x_);
	fftw_execute(planfor);
	for (int ii = 0; ii < nr; ++ii) {
		y[ii] = std::complex<double>(y_[ii][0], y_[ii][1]);
	}
	return;
}

/// The result gets stored in the object member `out` of type `double`, which is
/// a 1-d array.  Elements are accessed via e.g., `out[0]`, `out[n-1]`.
inline void VF::ifft(double* y, std::complex<double>* x) {
	for (int ii = 0; ii < nr; ++ii) {
		y_[ii][0] = real(x[ii]);		
		y_[ii][1] = imag(x[ii]);
	}

	fftw_execute(planback);
	for (int ii = 0; ii < n_; ++ii) {
		y[ii] = x_[ii] / n_;
	}
  return;
}

/// Product between `complex<double>` vectors.
///
/// Result is `z = x * y`.
///
/// @param[out] z Output complex vector, 2-d array.
/// @param[in] x First input complex vector, 2-d array.
/// @param[in] y Second input complex vector, 2-d array.
/// @param[in] n Size of inputs (integer).
inline void Complex_Mult(std::complex<double>* z, std::complex<double>* x, std::complex<double>* y,
                    int n) {
  for (int ii = 0; ii < n; ++ii) {
    z[ii] = x[ii] * y[ii];
  }
  return;
}

/// In-place subtract product between `complex<double>` vectors.
///
/// Result is `z -= x * y`.
///
/// @param[out] z Output complex vector, 2-d array.
/// @param[in] x First input complex vector, 2-d array.
/// @param[in] y Second input complex vector, 2-d array.
/// @param[in] n Size of inputs (integer).
inline void Complex_Mult_Minus(std::complex<double>* z, std::complex<double>* x,
	std::complex<double>* y, int n) {
  for (int ii = 0; ii < n; ++ii) {
	  z[ii] -= x[ii] * y[ii];
  }
  return;
}

/// In-place sum product between `complex<double>` vectors.
///
/// Result is `z += x * y`.
///
/// @param[out] z Output complex vector, 2-d array.
/// @param[in] x First input complex vector, 2-d array.
/// @param[in] y Second input complex vector, 2-d array.
/// @param[in] n Size of inputs (integer).
inline void Complex_Mult_Plus(std::complex<double>* z, std::complex<double>* x,
	std::complex<double>* y, int n) {
  for (int ii = 0; ii < n; ++ii) {
	  z[ii] += x[ii] * y[ii];
  }
  return;
}

/// Steps for polynomial convolution: For two polynomials `a(x) = a_0 + a_1 * x
/// + ... a_p * x^p` and  `b(x) = b_0 + b_1 * x + ... b_q * x^q`
///
/// 1) generate a VectorFFT classes, three vectors of complex<double>, the zero-padding vector of a,b and the storage space of c
///		v1 = new VectorFFT(p+q); 
///     std::complex<double>* a_fft = new complex<double>[p+q]
///     std::complex<double>* b_fft = new complex<double>[p+q]
///     std::complex<double>* c_fft = new complex<double>[p+q]
///     double* a0 = new double[p+q]
///     double* b0 = new double[p+q]
///     double* c = new double[p+q]
///
/// 2) copy the value of a and b into a0 and b0
///		std::copy(a, a + p, a0);
///		std::fill(a0+p, a0 + p + q, 0);
///		std::copy(b, b + q, b0);
///		std::fill(b0+q, b0 + p + q, 0);
///     
/// 3) compute the FFT of a0 and b0: 
///		v1->fft(a_fft, a0);
///		v1->fft(b_fft, b0);
///
/// 4) multiplication: 
///		Complex_Mult(c_fft, a_fft, b_fft, p+q);
///
/// 5) Inverse FFT: 
///		v1->ifft(c, c_fft);
///
/// Then results `c(x) = c_0 + c_1 * x + ... + c_p+q * x^p+q = a(x) * b(x)` is
/// in c.
///
/// Note 1: for two length-n polynomials a_n(x) and b_n(x), their VectorFFT is
/// of size 2*n. 
///
/// Note 2: In fftw, if c_FFT is of length 2*n and it is the FFT of
/// real vector, computation of its IFFT only requires its first 2*(N/2+1)
/// terms, not all of the 2*n terms.


void printVector(double* x, int n) {
	for (int ii = 0; ii < n; ++ii) {
		std::cout << x[ii] << " ";
	}
	std::cout << std::endl;
	return;
}

void printComplex(std::complex<double>* x, int n) {
	for (int ii = 0; ii < n; ++ii) {
		std::cout << x[ii] << " ";
	}
	std::cout << std::endl;
	return;
}


#endif

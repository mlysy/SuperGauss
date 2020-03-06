/// @file VectorFFT.h

///////////////////////////////////////////////
// Fast Fourier Transformation and Inverse
///////////////////////////////////////////////

#ifndef VectorFFT_h
#define VectorFFT_h 1

// usual header
#include <complex>
#include <fftw3.h>
// #include <Rcpp.h>
// #include <iostream>

/// FFT for real to complex and corresponding iFFT for complex to real.
///
/// @brief Allocates memory for the corresponding `fftw` operations within the object, copies memory in and out each time the FFT and iFFT members are called.
class VectorFFT {
private:
  typedef std::complex<double> complexd; ///< Typedef for complex double
  fftw_plan planfwd_;  ///< Plan for forward FFT.
  fftw_plan planback_;  ///< Plan for backward FFT.
  fftw_complex* y_; ///< Where to compute FFT.
  double* x_; ///< Where to compute iFFT.
  int n_; ///< Size of input vector.
  int npad_; ///< Padded size of input.
public:
  /// Constructor.
  VectorFFT(int n);
  /// Destructor.
  ~VectorFFT();
  /// Perform the FFT on the input data.
  void fft(std::complex<double>* y, const double* x);
  /// Perform the inverse FFT on the input data.
  void ifft(double* x, const std::complex<double>* y);
  /// Get size of FFT.
  int size() {return n_;}
};

/// @param[in] n Size of FFT/iFFT to be computed.
inline VectorFFT::VectorFFT(int n) {
  n_ = n;
  npad_ = ceil((double)(n + 1) / 2);
  x_ = new double[n];
  std::fill(x_, x_ + n, 0);
  y_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  planfwd_ = fftw_plan_dft_r2c_1d(n, x_, y_, FFTW_ESTIMATE);
  planback_ = fftw_plan_dft_c2r_1d(n, y_, x_, FFTW_ESTIMATE);
  return;
}

inline VectorFFT::~VectorFFT() {
  delete[] x_;
  fftw_free(y_);
  fftw_destroy_plan(planfwd_);
  fftw_destroy_plan(planback_);
}

/// Calculates `y = FFT(x)`.
///
/// @param[out] y Complex vector output.
/// @param[in] x  Real vector input.
inline void VectorFFT::fft(std::complex<double>* y,
			   const double* x) {
  std::copy(x, x + n_, x_);
  fftw_execute(planfwd_);
  for (int ii = 0; ii < npad_; ++ii) {
    y[ii] = complexd(y_[ii][0], y_[ii][1]);
  }
  return;
}

/// Calculates `x = iFFT(y)`.
///
/// @param[out] x Real vector output.
/// @param[in] y Complex vector input.
inline void VectorFFT::ifft(double* x, const std::complex<double>* y) {
  for (int ii = 0; ii < npad_; ++ii) {
    y_[ii][0] = real(y[ii]);
    y_[ii][1] = imag(y[ii]);
  }

  fftw_execute(planback_);
  for (int ii = 0; ii < n_; ++ii) {
    x[ii] = x_[ii] / n_;
  }
  return;
}

/// Elementwise product between `std::complex<double>` arrays.
///
/// Result is `z[i] = x[i] * y[i]`.
///
/// @param[out] z Output complex array.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of each array.
inline void Complex_Mult(std::complex<double>* z,
			 const std::complex<double>* x,
			 const std::complex<double>* y,
			 int n) {
  for (int ii = 0; ii < n; ++ii) {
    z[ii] = x[ii] * y[ii];
  }
  return;
}

/// Elementwise subtract-product between `std::complex<double>` vectors.
///
/// Result is `z[i] -= x[i] * y[i]`.
///
/// @param[out] z Output complex array.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of each array.
inline void Complex_Mult_Minus(std::complex<double>* z,
			       const std::complex<double>* x,
			       const std::complex<double>* y,
			       int n) {
  for (int ii = 0; ii < n; ++ii) {
    z[ii] -= x[ii] * y[ii];
  }
  return;
}

/// Elementwise add-product between `std::complex<double>` arrays.
///
/// Result is `z[i] += x[i] * y[i]`.
///
/// @param[out] z Output complex array.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of each array.
inline void Complex_Mult_Plus(std::complex<double>* z,
			      const std::complex<double>* x,
			      const std::complex<double>* y,
			      int n) {
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


#endif

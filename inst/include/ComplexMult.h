/// @file ComplexMult.h
/// @brief Helper functions for polynomial convolutions.
///
/// Consider two polynomials `a(x) = a_0 + a_1 * x + ... a_p * x^p` and `b(x) = b_0 + b_1 * x + ... b_q * x^q`, of which the coefficients are stored as `double` arrays `a` and `b`.  Now suppose we wish to compute the coefficients of `c(x) = a(x) * b(x)` and store them in another `double` array `c`.  Then we have the following steps: 
///
/// 1.  Generate a `VectorFFT` object, three arrays of `std:complex<double>`, the zero-padding vector of `a`, and `b` and the storage space of for `c`:
///
///     ```
///     v1 = new VectorFFT(p+q); 
///     std::complex<double>* a_fft = new std::complex<double>[p+q]
///     std::complex<double>* b_fft = new std::complex<double>[p+q]
///     std::complex<double>* c_fft = new std::complex<double>[p+q]
///     double* a0 = new double[p+q]
///     double* b0 = new double[p+q]
///     double* c = new double[p+q]
///     ```
///
/// 2. Copy the value of `a` and `b` into `a0` and `b0`:
///
///    ```
///    std::copy(a, a + p, a0);
///    std::fill(a0+p, a0 + p + q, 0);
///    std::copy(b, b + q, b0);
///    std::fill(b0+q, b0 + p + q, 0);
///    ```
///     
/// 3. Compute the FFT of `a0` and `b0`:
///
///     ```
///     v1->fft(a_fft, a0);
///     v1->fft(b_fft, b0);
///     ```
///
/// 4. Do elementwise multiplication in the Fourier domain:
///
///    ```
///    Complex_Mult(c_fft, a_fft, b_fft, p+q);
///    ```
///
/// 5. Inverse FFT to recover `c`:
///
///     ```
///     v1->ifft(c, c_fft);
///     ```
///
/// @note
/// 1. For two length-`n` polynomials `a_n(x)` and `b_n(x)`, their `VectorFFT` is of size `2*n`. 
///
/// 2. In `fftw`, if `c_fft` is of length `2*n` and it is the FFT of a real vector, computation of its iFFT only requires its first `2*(n/2+1)` terms, not all `2*n` terms.

#ifndef ComplexMult_h
#define ComplexMult_h 1


/// Elementwise product between `std::complex<double>` arrays.
///
/// Result is `z[i] = x[i] * y[i]`.
///
/// @param[out] z Output complex array.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of each array.
inline void complex_mult(std::complex<double>* z,
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
inline void complex_mult_minus(std::complex<double>* z,
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
inline void complex_mult_plus(std::complex<double>* z,
			      const std::complex<double>* x,
			      const std::complex<double>* y,
			      int n) {
  for (int ii = 0; ii < n; ++ii) {
    z[ii] += x[ii] * y[ii];
  }
  return;
}

#endif

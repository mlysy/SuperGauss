/// @file VectorFFT.h

///////////////////////////////////////////////
// Fast Fourier Transformation and Inverse
///////////////////////////////////////////////

#ifndef VectorFFT_h
#define VectorFFT_h 1

// usual header
#include <Rcpp.h>
#include <fftw3.h>
#include <iostream>
#include <ctime>
// using namespace Rcpp;
// using namespace std;

/// Forward FFT from real to complex.
class VectorFFT {
 private:
  fftw_plan planfor;  ///< FFTW plan.
 public:
  int n_size;         ///< Size of input vector.
  double* in;         ///< Real input vector.
  fftw_complex* out;  ///< Complex output vector.

  /// Constructor.
  VectorFFT(int);
  /// Destructor.
  ~VectorFFT();
  /// Perform the FFT on the input data.
  void fft();
};

/// @param[in] n Size of the input vector.
inline VectorFFT::VectorFFT(int n) {
  n_size = n;
  in = new double[n];
  out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  std::fill(in, in + n, 0);
  planfor = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
  return;
}

inline VectorFFT::~VectorFFT() {
  delete[] in;
  fftw_free(out);
  fftw_destroy_plan(planfor);
}

/// The result gets stored in the object member `out` of type `fftw_complex`,
/// which is a 2-d array.  Elements are accessed via e.g., `out[0][1]`,
/// `out[n-1][0]`.
inline void VectorFFT::fft() {
  fftw_execute(planfor);
  return;
}

/// Backward FFT from complex to real.
class VectorIFFT {
  fftw_plan planback;  ///< FFTW plan.
 public:
  int n_size;        ///< Size of input vector.
  double* out;       ///< Real output vector.
  fftw_complex* in;  ///< Complex input vector.
  /// Constructor.
  VectorIFFT(int);
  /// Perform the inverse FFT on the input data.
  void Ifft();
  /// Destructor.
  ~VectorIFFT();
};

/// @param[in] n Size of the input vector.
inline VectorIFFT::VectorIFFT(int n) {
  n_size = n;
  out = new double[n];
  in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  memset(in, 0, sizeof(fftw_complex) * n);
  planback = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
}

inline VectorIFFT::~VectorIFFT() {
  delete[] out;
  fftw_free(in);
  fftw_destroy_plan(planback);
}

/// The result gets stored in the object member `out` of type `double`, which is
/// a 1-d array.  Elements are accessed via e.g., `out[0]`, `out[n-1]`.
inline void VectorIFFT::Ifft() {
  fftw_execute(planback);
  for (int ii = 0; ii < n_size; ++ii) {
    out[ii] = out[ii] / n_size;
  }
  return;
}

/// Product between `fftw_complex` vectors.
///
/// Result is `y = alpha * beta`.
///
/// @param[out] y Output complex vector, 2-d array.
/// @param[in] alpha First input complex vector, 2-d array.
/// @param[in] beta Second input complex vector, 2-d array.
/// @param[in] n Size of inputs (integer).
inline void vecConv(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta,
                    int n) {
  for (int ii = 0; ii < n; ++ii) {
    y[ii][0] = alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
    y[ii][1] = alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
  }
  return;
}

/// In-place subtract product between `fftw_complex` vectors.
///
/// Result is `y -= alpha * beta`.
///
/// @param[out] y Output complex vector, 2-d array.
/// @param[in] alpha First input complex vector, 2-d array.
/// @param[in] beta Second input complex vector, 2-d array.
/// @param[in] n Size of inputs (integer).
inline void vecConv_Sub(fftw_complex* y, fftw_complex* alpha,
                        fftw_complex* beta, int n) {
  for (int ii = 0; ii < n; ++ii) {
    y[ii][0] -= alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
    y[ii][1] -= alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
  }
  return;
}

/// In-place sum product between `fftw_complex` vectors.
///
/// Result is `y += alpha * beta`.
///
/// @param[out] y Output complex vector, 2-d array.
/// @param[in] alpha First input complex vector, 2-d array.
/// @param[in] beta Second input complex vector, 2-d array.
/// @param[in] n Size of inputs (integer).
inline void vecConv_Add(fftw_complex* y, fftw_complex* alpha,
                        fftw_complex* beta, int n) {
  for (int ii = 0; ii < n; ++ii) {
    y[ii][0] += alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
    y[ii][1] += alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
  }
  return;
}

/// Steps for polynomial convolution: For two polynomials `a(x) = a_0 + a_1 * x
/// + ... a_p * x^p` and  `b(x) = b_0 + b_1 * x + ... b_q * x^q`
///
/// 1) generate two VectorFFT classes: a_FFT = new VectorFFT(p+q), b_FFT = new
/// VectorFFT(p+q) and one VectorIFFT class: c_IFFT = new VectorIFFT(p+q)
///
/// 2) fill the VectorFFT->in: std::copy(a, a+p, a_FFT->in), std::copy(b, b+q,
/// b_FFT->in)
///
/// 3) FFT: a_FFT->fft(), b_FFT->fft()
///
/// 4) multiplication between VectorFFT->out: vecConv(c_FFT->out, a_FFT->out,
/// b_FFT->out, p+q)
///
/// 5) IFFT: c_IFFT->ifft()
///
/// Then results `c(x) = c_0 + c_1 * x + ... + c_p+q * x^p+q = a(x) * b(x)` is
/// in c_IFFT->out.
///
/// Note 1: for two length-n polynomials a_n(x) and b_n(x), their VectorFFT is
/// of size 2*n. 
///
/// Note 2: In fftw, if c_FFT is of length 2*n and it is the FFT of
/// real vector, computation of its IFFT only requires its first 2*(N/2+1)
/// terms, not all of the 2*n terms.

#endif

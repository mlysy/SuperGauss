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
// #include <ctime> 
using namespace Rcpp; // FIXME: Never put "using namespace" in header file
using namespace std; // REMOVE!

/// Notes: Usage of VectorFFT class
///
/// For two polynomials a^(p)(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_p * x^p and b^(q)(x) = b_0 + b_1 * x + b_2 * x^2 + ... b_q * x^q
///
/// Their convolution result c^(p+q)(x) = c_0 + c_1 * x + c_2 * x^2 + ... + c_{p+q} * x^{p+q} = conv(a(x), b(x)) can be computed efficiently using Fast Fourier Transformation in following steps:
///
/// 1) Define two length-(p+q) vectors: aa = (a_0, a_1, ..., a_p, 0, ..., 0) and bb = (b_0, b_1, ..., b_q, 0, ..., 0).
///
/// 2) Compute their forward-fft: aa_fft = fft(aa) and bb_fft = fft(b), both are 2-d array of complex vector. First array aa_fft[:][0] and bb_fft[:][0] are the real parts and the second array array aa_fft[:][1] and bb_fft[:][1] are the complex parts.
///
/// 3) Multiply: cc_fft = aa_fft * bb_fft, where cc_fft is also a 2-d array of complex vector. In c++ it is computed like this: cc_fft[:][0] = aa_fft[:][0] * bb_fft[:][0] - aa_fft[:][1] * bb_fft[:][1], cc_fft[:][1] = aa_fft[:][0] * bb_fft[:][1] + aa_fft[:][1] * bb_fft[:][0]
///
/// 4) Backward-fft: ifft(cc_fft) to obtain the length p+q vector cc = (c_0, c_1, ..., c_{p+q})

/// Forward FFT from real to complex.
class VectorFFT{
private:
  fftw_plan planfor; ///< FFTW plan.
public:
  int n_size; ///< Size of input vector.
  double* in; ///< Real input vector.
  fftw_complex *out; ///< Complex output vector.
  /// Constructor.
  VectorFFT(int);
  /// Perform the FFT on the input data.
  void fft();
  /// Destructor.
  ~VectorFFT();
};

/// @param[in] n Size of the input vector.
inline VectorFFT::VectorFFT(int n){
  n_size = n;
  in = new double[n];
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  std::fill(in, in + n, 0);
  planfor = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
  return;
}

inline VectorFFT::~VectorFFT(){
  delete []in;
  fftw_free(out);
  fftw_destroy_plan(planfor);
}

/// The result gets stored in the object member `out` of type `fftw_complex`, which is a 2-d array.  Elements are accessed via e.g., `out[0][1]`, `out[n-1][0]`.
inline void VectorFFT::fft(){
  fftw_execute(planfor);
  return;
}

/// Backward FFT from complex to real.
class VectorIFFT{
  fftw_plan planback; ///< FFTW plan.
public:
  int n_size; ///< Size of input vector.
  double* out; ///< Real output vector.
  fftw_complex *in; ///< Complex input vector.
  /// Constructor.
  VectorIFFT(int);
  /// Perform the inverse FFT on the input data.
  void Ifft();
  /// Destructor.
  ~VectorIFFT();
};

/// @param[in] n Size of the input vector.
inline VectorIFFT::VectorIFFT(int n){
  n_size = n;
  out = new double[n];
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  memset(in, 0, sizeof(fftw_complex) * n);
  planback = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
}

inline VectorIFFT::~VectorIFFT(){
  delete []out;
  fftw_free(in);
  fftw_destroy_plan(planback);
}

/// The result gets stored in the object member `out` of type `double`, which is a 1-d array.  Elements are accessed via e.g., `out[0]`, `out[n-1]`.
inline void VectorIFFT::Ifft(){
  fftw_execute(planback);
  for(int ii = 0;ii < n_size; ++ii){
    out[ii] = out[ii]/n_size;
  }
  return;
}

/// Product between `fftw_complex` vectors.
///
/// Result is `y = alpha * beta`.
///
/// @param[out] y Output vector.
/// @param[in] alpha First input vector.
/// @param[in] beta Second input vector.
/// @param[in] n Size of inputs (integer).
inline void vecConv(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n){
	for (int ii = 0; ii < n; ++ii){
		y[ii][0] = alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] = alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

/// In-place subtract product between `fftw_complex` vectors.
///
/// Result is `y -= alpha * beta`.
///
/// @param[out] y Output vector.
/// @param[in] alpha First input vector.
/// @param[in] beta Second input vector.
/// @param[in] n Size of inputs (integer).
inline void vecConv_Sub(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n){
	for (int ii = 0; ii < n; ++ii){
		y[ii][0] -= alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] -= alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

/// In-place sum product between `fftw_complex` vectors.
///
/// Result is `y += alpha * beta`.
///
/// @param[out] y Output vector.
/// @param[in] alpha First input vector.
/// @param[in] beta Second input vector.
/// @param[in] n Size of inputs (integer).
inline void vecConv_Add(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n){
	for (int ii = 0; ii < n; ++ii){
		y[ii][0] += alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] += alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

# endif

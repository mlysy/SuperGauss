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

///////////////////////////////////////////////
// Fast Fourier Transformation and Inverse
///////////////////////////////////////////////

#ifndef VectorFFT_h
#define VectorFFT_h 1

// usual header
#include <Rcpp.h>
using namespace Rcpp;
#include <fftw3.h>
#include <iostream>
#include <ctime>
using namespace std;


// defining classes
//------------------------------------------------------
// 1, fast fourier transformation
class VectorFFT{
  fftw_plan planback;
public:
  int n_size;
  double* in;
  fftw_complex *out;
  VectorFFT(int);
  void fft();
  ~VectorFFT();
};

//------------------------------------------------------
// 1, fast fourier transformation
// class VectorFFT

inline VectorFFT::VectorFFT(int n){
  n_size = n;
  in = new double[n];
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  std::fill(in, in + n, 0);
  planback = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
  return;
}

inline VectorFFT::~VectorFFT(){
  delete []in;
  fftw_free(out);
  fftw_destroy_plan(planback);
}

inline void VectorFFT::fft(){
  fftw_execute(planback);
  return;
}
//------------------------------------------------------

//------------------------------------------------------

// 2, inverse fast fourier transformation
class VectorIFFT{
  fftw_plan planfor;
public:
  int n_size;
  double* out;
  fftw_complex *in;
  VectorIFFT(int);
  void Ifft();
  ~VectorIFFT();
};
//------------------------------------------------------

// 2, inverse fast fourier transformation
//  class VectorIFFT

inline VectorIFFT::VectorIFFT(int n){
  n_size = n;
  out = new double[n];
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  memset(in, 0, sizeof(fftw_complex) * n);
  planfor = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
}

inline VectorIFFT::~VectorIFFT(){
  delete []out;
  fftw_free(in);
  fftw_destroy_plan(planfor);
}

inline void VectorIFFT::Ifft(){
  fftw_execute(planfor);
  for(int ii = 0;ii < n_size; ++ii){
    out[ii] = out[ii]/n_size;
  }
  return;
}
//------------------------------------------------------

inline void vecConv(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n){
	for (int ii = 0; ii < n; ++ii){
		y[ii][0] = alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] = alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

inline void vecConv_Sub(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n){
	for (int ii = 0; ii < n; ++ii){
		y[ii][0] -= alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] -= alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

inline void vecConv_Add(fftw_complex* y, fftw_complex* alpha, fftw_complex* beta, int n){
	for (int ii = 0; ii < n; ++ii){
		y[ii][0] += alpha[ii][0] * beta[ii][0] - alpha[ii][1] * beta[ii][1];
		y[ii][1] += alpha[ii][1] * beta[ii][0] + alpha[ii][0] * beta[ii][1];
	}
	return;
}

# endif
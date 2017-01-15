///////////////////////////////////////////////
// Fast Fourier Transformation and Inverse
///////////////////////////////////////////////

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
class VectorFFT 
{
  fftw_plan planback;
public:
  double* in;
  fftw_complex *out;
  VectorFFT(int);
  void fft();
  ~VectorFFT();
};
//------------------------------------------------------

// 2, inverse fast fourier transformation
class VectorIFFT 
{
  fftw_plan planfor;
  int n_size;
public:
  double* out;
  fftw_complex *in;
  VectorIFFT(int);
  void Ifft();
  ~VectorIFFT();
};
//------------------------------------------------------
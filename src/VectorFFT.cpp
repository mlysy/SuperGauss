//defining functions in different classes

#include "VectorFFT.h"

//------------------------------------------------------
// 1, fast fourier transformation
// class VectorFFT

VectorFFT::VectorFFT(int n)
{
  in = new double[n];
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  std::fill(in, in + n, 0);
  planback = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
  return;
}

VectorFFT::~VectorFFT()
{
  delete []in;
  fftw_free(out);
  fftw_destroy_plan(planback);
}

void VectorFFT::fft()
{
  fftw_execute(planback);
  return;
}
//------------------------------------------------------

// 2, inverse fast fourier transformation
//  class VectorIFFT

VectorIFFT::VectorIFFT(int n)
{
  n_size = n;
  out = new double[n];
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  memset(in, 0, sizeof(fftw_complex) * n);
  planfor = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
}

VectorIFFT::~VectorIFFT()
{
  delete []out;
  fftw_free(in);
  fftw_destroy_plan(planfor);
}

void VectorIFFT::Ifft()
{
  fftw_execute(planfor);
  for(int ii = 0;ii < n_size; ++ii)
  {
    out[ii] = out[ii]/n_size;
  }
  return;
}
//------------------------------------------------------
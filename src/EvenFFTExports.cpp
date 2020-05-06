/// @file EvenFFTExports.cpp

#include "EvenFFT.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector even_fft(NumericVector x, bool inverse = false) {
  int N = x.length();
  EvenFFT efft(N);
  int Nu = efft.usize(); // storage size
  NumericVector y(Nu);
  if(!inverse) {
    efft.fft(REAL(y), REAL(x));
  } else {
    efft.ifft(REAL(y), REAL(x));
  }
  return y;
}

// enum EFFT_type {E00, E01, E10, E11, O00, O01, O10, O11};

// class EFFT {
// private:
//   int n_;
//   double *x_;
//   double *y_;
//   fftw_plan plan_;
// public:
//   EFFT(int n, EFFT_type type) {
//     n_ = n;
//     x_ = fftw_alloc_real(n_);
//     y_ = fftw_alloc_real(n_);
//     switch(type) {
//       case E00 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_REDFT00, FFTW_ESTIMATE);
//       case E01 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_REDFT01, FFTW_ESTIMATE);
//       case E10 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_REDFT10, FFTW_ESTIMATE);
//       case E11 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_REDFT11, FFTW_ESTIMATE);
//       case O00 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_RODFT00, FFTW_ESTIMATE);
//       case O01 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_RODFT01, FFTW_ESTIMATE);
//       case O10 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_RODFT10, FFTW_ESTIMATE);
//       case O11 :
// 	plan_ = fftw_plan_r2r_1d(n_, x_, y_, FFTW_RODFT11, FFTW_ESTIMATE);
//     }
//   }
//   ~EFFT() {
//     fftw_free(x_);
//     fftw_free(y_);
//     fftw_destroy_plan(plan_);
//   }
//   void fwd(double* y, const double* x) {
//     std::copy(x, x+n_, x_);
//     fftw_execute(plan_);
//     std::copy(y_, y_+n_, y);
//     return;
//   }
// };

// // [[Rcpp::export]]
// NumericVector even_fft(NumericVector x, int type) {
//   int n = x.length();
//   NumericVector y(n);
//   EFFT* efft;
//   switch(type) {
//     case 0 :
//       efft = new EFFT(n, E00);
//     case 1 :
//       efft = new EFFT(n, E01);
//     case 2 :
//       efft = new EFFT(n, E10);
//     case 3 :
//       efft = new EFFT(n, E11);
//     case 4 :
//       efft = new EFFT(n, O00);
//     case 5 :
//       efft = new EFFT(n, O01);
//     case 6 :
//       efft = new EFFT(n, O10);
//     case 7 :
//       efft = new EFFT(n, O11);
//   }
//   efft->fwd(REAL(y), REAL(x));
//   delete efft;
//   return y;
// }

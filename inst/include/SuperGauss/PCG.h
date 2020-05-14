/// @file PCG.h

#ifndef PCG_h
#define PCG_h 1

#include "RealFFT.h"

double crossprod(double* a, double* b, int n) {
  double pd = 0;
  for (int ii = 0; ii < n; ++ii) {
    pd += a[ii] * b[ii];
  }
  return pd;
}

class PCG {
private:
  int N; // size
  double* pchan; // storage of pre conditional
  RealFFT* tf; // Toeplitz fft, size 2*N
  RealFFT* cf; // circulant fft, size N
  double* x;
  std::complex<double>* fx;
  double* x2;
  std::complex<double>* fx2;
  std::complex<double>* facf;
  std::complex<double>* ifcv;
  double* tm; // store the result of toep_mult
  double* rr;
  double* zz; // store the result of circ_solve
  double* dd;

public:
  // Constructor
  PCG(int);
  // Destructor
  ~PCG();
  // Toeplitz mult
  void toep_mult(double*);
  // Circulant solve
  void circ_solve(double*);
  // Solve
  void solve(double*, double*, double*, double);
};

inline PCG::PCG(int n) {
  N = n;
  pchan = new double[n];
  tf = new RealFFT(2 * n);
  cf = new RealFFT(n);
  x = new double[2 * N];
  fx = new std::complex<double>[N];
  x2 = new double[2 * N];
  fx2 = new std::complex<double>[2 * N];
  facf = new std::complex<double>[2 * N];
  ifcv = new std::complex<double>[N];
  tm = new double[n];
  rr = new double[n];
  zz = new double[n];
  dd = new double[n];
}

inline PCG::~PCG() {
  delete[] pchan;
  delete tf;
  delete cf;
  delete[] x;
  delete[] fx;
  delete[] x2;
  delete[] fx2;
  delete[] facf;
  delete[] ifcv;
  delete[] tm;
  delete[] rr;
  delete[] zz;
  delete[] dd;
}

inline void PCG::toep_mult(double* in) {
  std::copy(in, in + N, x2);
  std::fill(x2 + N, x2 + 2 * N, 0);
  tf->fft(fx2, x2);
  for (int ii = 0; ii < 2 * N; ++ii) {
    fx2[ii] *= facf[ii];
  }
  tf->ifft(x2, fx2);
  std::copy(x2, x2 + N, tm);
}


inline void PCG::circ_solve(double* in) {
  std::copy(in, in + N, x);
  cf->fft(fx, x);
  for (int ii = 0; ii < N; ++ii) {
    fx[ii] *= ifcv[ii];
  }
  cf->ifft(zz, fx);
}

inline void PCG::solve(double* xx, double* acf, double* y, double tol) {
  double cst = 2 * 3.14159265358979323846 * N;
  pchan[0] = N * acf[0] / cst;
  for (int ii = 1; ii < N; ++ii) {
    pchan[ii] = (N - ii) * acf[ii] + ii * acf[N - ii];
    pchan[ii] /= cst;
  }
  std::copy(acf, acf + N, x);
  x[N] = 0;
  std::copy(acf + 1, acf + N, x + N + 1);
  std::reverse(x + N + 1, x + 2 * N);
  tf->fft(facf, x);
  cf->fft(ifcv, pchan);
  for (int ii = 0; ii < N; ++ii) {
    ifcv[ii] = 1. / ifcv[ii];
  }
  std::fill(xx, xx + N, 0);
  std::copy(y, y + N, rr);
  circ_solve(rr);
  std::copy(zz, zz + N, dd);
  double rz = crossprod(rr, zz, N);
  double alpha, beta;
  for (int ii = 0; ii < N; ++ii) {
    if (sqrt(crossprod(rr, rr, N)) < tol) {
      break;
    }
    toep_mult(dd);
    alpha = rz / crossprod(dd, tm, N);
    for (int jj = 0; jj < N; jj++) {
      xx[jj] += alpha * dd[jj];
      rr[jj] -= alpha * tm[jj];
    }
    circ_solve(rr);
    beta = rz;
    rz = crossprod(rr, zz, N);
    beta = rz / beta;
    for (int jj = 0; jj < N; ++jj) {
      dd[jj] = zz[jj] + beta * dd[jj];
    }
  }
}

#endif

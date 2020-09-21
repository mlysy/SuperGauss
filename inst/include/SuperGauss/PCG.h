/// @file PCG.h

#ifndef PCG_h
#define PCG_h 1

#include "RealFFT.h"

class PCG {
private:
  typedef std::complex<double> dcomplex;
  int N_; // size
  double* pchan_; // storage of pre conditional
  RealFFT* tfft_; // Toeplitz fft, size 2*N_
  RealFFT* cfft_; // circulant fft, size N_
  double* res_mult_; // store the result of toep_mult
  dcomplex* facf_; // used in toep_mult
  double* res_solve_; // store the result of circ_solve
  dcomplex* ifcv_; // used in circ_solve
  double* vec1_;
  dcomplex* vec1_fft_;
  double* vec2_;
  dcomplex* vec2_fft_;
  double* rr_pcg_;
  double* dd_pcg_;
  /// Dot-product between vectors.
  double dot_prod(const double* v1, const double* v2);
  /// Internal Toeplitz matrix multiplication.
  void toep_mult(const double* x);
  /// Internal Circulant matrix linear system solver.
  void circ_solve(const double* x);

public:
  /// Constructor.
  PCG(int N);
  /// Destructor.
  ~PCG();
  /// Solve Toeplitz matrix-vector system of equations.
  void solve(double* y, const double* acf, const double* x, double tol);
};

/// @param[in] N Size of the Toeplitz matrix.
inline PCG::PCG(int N) {
  N_ = N;
  pchan_ = new double[N_];
  tfft_ = new RealFFT(2 * N_);
  cfft_ = new RealFFT(N_);
  vec1_ = new double[2 * N_];
  vec1_fft_ = new dcomplex[N_];
  vec2_ = new double[2 * N_];
  vec2_fft_ = new dcomplex[2 * N_];
  facf_ = new dcomplex[2 * N_];
  ifcv_ = new dcomplex[N_];
  res_mult_ = new double[N_];
  res_solve_ = new double[N_];
  rr_pcg_ = new double[N_];
  dd_pcg_ = new double[N_];
}

inline PCG::~PCG() {
  delete[] pchan_;
  delete tfft_;
  delete cfft_;
  delete[] vec1_;
  delete[] vec1_fft_;
  delete[] vec2_;
  delete[] vec2_fft_;
  delete[] facf_;
  delete[] ifcv_;
  delete[] res_mult_;
  delete[] res_solve_;
  delete[] rr_pcg_;
  delete[] dd_pcg_;
}

/// @param[in] v1 First vector of length `N`.
/// @param[in] v2 Second vector of length `N`.
/// @return Dot product `sum_i=1^n v1_i * v2_i`.
inline double PCG::dot_prod(const double* v1, const double* v2) {
  double ans = 0.0;
  for(int ii=0; ii<N_; ++ii) {
    ans += v1[ii] * v2[ii];
  }
  return ans;
}

/// @param[in] x Input vector of length `N`.
inline void PCG::toep_mult(const double* x) {
  // temporaries: vec2_, vec2_fft_
  std::copy(x, x + N_, vec2_);
  std::fill(vec2_ + N_, vec2_ + 2 * N_, 0);
  tfft_->fft(vec2_fft_, vec2_);
  for(int ii = 0; ii < 2 * N_; ++ii) {
    vec2_fft_[ii] *= facf_[ii];
  }
  tfft_->ifft(vec2_, vec2_fft_);
  std::copy(vec2_, vec2_ + N_, res_mult_);
  return;
}

/// @param[in] x Input vector of length `N`.
inline void PCG::circ_solve(const double* x) {
  // temporaries: vec1_, vec1_fft_
  std::copy(x, x + N_, vec1_);
  cfft_->fft(vec1_fft_, vec1_);
  for(int ii = 0; ii < N_; ++ii) {
    vec1_fft_[ii] *= ifcv_[ii];
  }
  cfft_->ifft(res_solve_, vec1_fft_);
  return;
}

/// @param[in] y Vector of length `N` in which to store the solution.
/// @param[in] acf Vector of length `N` containing the autocorrelation.
/// @param[in] x Vector of length `N` containing the input to the system.
/// @param[in] tol Nonnegative scalar specifying the tolerance for the conjugate gradient.
inline void PCG::solve(double* y, const double* acf, const double* x,
		       double tol) {
  // temporaries: pchan_, vec1_, rr_pcg_, dd_pcg_
  double cst = 2 * 3.14159265358979323846 * N_;
  double rz, alpha, beta;
  pchan_[0] = N_ * acf[0] / cst;
  for(int ii = 1; ii < N_; ++ii) {
    pchan_[ii] = (N_ - ii) * acf[ii] + ii * acf[N_ - ii];
    pchan_[ii] /= cst;
  }
  std::copy(acf, acf + N_, vec1_);
  vec1_[N_] = 0;
  // std::copy(acf + 1, acf + N_, vec1_ + N_ + 1);
  // std::reverse(vec1_ + N_ + 1, vec1_ + 2 * N_);
  std::reverse_copy(acf + 1, acf + N_, vec1_ + N_ + 1);
  tfft_->fft(facf_, vec1_);
  cfft_->fft(ifcv_, pchan_);
  for(int ii = 0; ii < N_; ++ii) {
    ifcv_[ii] = 1. / ifcv_[ii];
  }
  std::fill(y, y + N_, 0);
  std::copy(x, x + N_, rr_pcg_);
  circ_solve(rr_pcg_);
  std::copy(res_solve_, res_solve_ + N_, dd_pcg_);
  rz = dot_prod(rr_pcg_, res_solve_);
  for(int ii = 0; ii < N_; ++ii) {
    if(sqrt(dot_prod(rr_pcg_, rr_pcg_)) < tol) {
      break;
    }
    toep_mult(dd_pcg_);
    alpha = rz / dot_prod(dd_pcg_, res_mult_);
    for(int jj = 0; jj < N_; jj++) {
      y[jj] += alpha * dd_pcg_[jj];
      rr_pcg_[jj] -= alpha * res_mult_[jj];
    }
    circ_solve(rr_pcg_);
    beta = rz;
    rz = dot_prod(rr_pcg_, res_solve_);
    beta = rz / beta;
    for(int jj = 0; jj < N_; ++jj) {
      dd_pcg_[jj] = res_solve_[jj] + beta * dd_pcg_[jj];
    }
  }
  return;
}

#endif

/// @file GSchur.h

///////////////////////////////////////////////
// Generalized Schur Algorithm
///////////////////////////////////////////////

#ifndef GSchur_h
#define GSchur_h 1

#include "RealFFT.h"
// #include "ComplexMult.h"
#include <vector>
#include <complex>
#include <algorithm>

/// @brief Memory allocation for the Generalized Schur Algorithm.
///
/// The GSchur algorithm involves several polynomial multiplications, for which the `GSchur2K` class provides storage to compute efficiently via the FFT.  The notation follows that of Ling & Lysy (2020).
///
/// @note Size of input is `2*n` for original GSchur algorithm and `n+m` for the merge algorithm.
struct GSchur2K {
private:
  typedef std::complex<double> dcomplex;
public:
  RealFFT* FFT; ///< Memory for the FFT itself.
  double* alpha0; ///< Coefficients of polynomial \f$\alpha_{0,2n}(x)\f$.
  dcomplex* alpha0_fft; ///< FFT of `alpha0`.
  double* alphan; ///< Coefficients of polynomial \f$\alpha_{n,n}(x)\f$.
  dcomplex* alphan_fft; ///< FFT of `alphan`.
  double* beta0; ///< Coefficients of polynomial \f$\beta_{0,2n}(x)\f$.
  dcomplex* beta0_fft; ///< FFT of `beta0`.
  double* betan; ///< Coefficients of polynomial \f$\beta_{n,n}(x)\f$.
  dcomplex* betan_fft; ///< FFT of `betan`.
  double* eta0; ///< Coefficients of polynomial \f$\eta_{0,n}(x)\f$.
  dcomplex* eta0_fft; ///< FFT of `eta0`.
  double* etat; ///< Coefficients of polynomial \f$\tilde\eta_{0,n}(x)\f$.
  dcomplex* etat_fft; ///< FFT of `etat`;
  double* etan; ///< Coefficients of polynomial \f$\eta_{n,n}(x)\f$.
  dcomplex* etan_fft; ///< FFT of `etan`.
  double* eta2n; ///< Coefficients of polynomial \f$\eta_{0,2n}(x)\f$.
  dcomplex* eta2n_fft; ///< FFT of `eta2n`.
  double* xi0; ///< Coefficients of polynomial \f$\xi_{0,n}(x)\f$.
  dcomplex* xi0_fft; ///< FFT of `xi0`.
  double* xit; ///< Coefficients of polynomial \f$\tilde\xi_{0,n}(x)\f$.
  dcomplex* xit_fft; ///< FFT of `xit`.
  double* xin; ///< Coefficients of polynomial \f$\xi_{n,n}(x)\f$.
  dcomplex* xin_fft; ///< FFT of `xin`.
  double* xi2n; ///< Coefficients of polynomial \f$\xi_{0,2n}(x)\f$.
  dcomplex* xi2n_fft; ///< FFT of `xi2n`.
  double* gamma;  ///< Coefficients of polynomial \f$\gamma(x)\f$.
  /// Constructor.
  GSchur2K(int n);
  /// Destructor.
  ~GSchur2K();
};

/// @param[in] n Size of each GSchur polynomial.
inline GSchur2K::GSchur2K(int n) {
  FFT = new RealFFT(n);
  alpha0 = new double[n];
  std::fill(alpha0, alpha0 + n, 0);
  alpha0_fft = new dcomplex[n];
  alphan = new double[n];
  alphan_fft = new dcomplex[n];
  beta0 = new double[n];
  std::fill(beta0, beta0 + n, 0);
  beta0_fft = new dcomplex[n];
  betan = new double[n];
  betan_fft = new dcomplex[n];
  eta0 = new double[n];
  std::fill(eta0, eta0 + n, 0);
  eta0_fft = new dcomplex[n];
  etat = new double[n];
  std::fill(etat, etat + n, 0);
  etat_fft = new dcomplex[n];
  etan = new double[n];
  etan_fft = new dcomplex[n];
  eta2n = new double[n];
  eta2n_fft = new dcomplex[n];
  xi0 = new double[n];
  std::fill(xi0, xi0 + n, 0);
  xi0_fft = new dcomplex[n];
  xit = new double[n];
  std::fill(xit, xit + n, 0);
  xit_fft = new dcomplex[n];
  xin = new double[n];
  xin_fft = new dcomplex[n];
  xi2n = new double[n];
  xi2n_fft = new dcomplex[n];
  gamma = new double[n];
}

inline GSchur2K::~GSchur2K() {
  delete FFT;
  delete[] alpha0;
  delete[] alpha0_fft;
  delete[] alphan;
  delete[] alphan_fft;
  delete[] beta0;
  delete[] beta0_fft;
  delete[] betan;
  delete[] betan_fft;
  delete[] eta0;
  delete[] eta0_fft;
  delete[] etat;
  delete[] etat_fft;
  delete[] etan;
  delete[] etan_fft;
  delete[] eta2n;
  delete[] eta2n_fft;
  delete[] xi0;
  delete[] xi0_fft;
  delete[] xit;
  delete[] xit_fft;
  delete[] xin;
  delete[] xin_fft;
  delete[] xi2n;
  delete[] xi2n_fft;
  delete[] gamma;
}

/// @brief The Generalized Schur Algorithm.
///
/// Given a symmetric positive-definite `N x N` Toeplitz matrix `Tz = Toeplitz(acf)` with first row/column `acf`, calculates `delta`, the first row/column of `Tz^{-1}` and `log(det(Tz))`.  This is done using a modified version of the "superfast" `O(N log^2(N))` algorithm of Ammar & Gragg (1988) which consists of the following steps:
///
/// 1. Decompose `N` into a vectors `s = (s_1, ..., s_k)`, where each `s_i` is a power of 2 and `sum(s) = N`.
///
/// 2. Divide the problem into pieces of size `s_i`, and solve each using the recursive algorithm of Ammar & Gragg (1988).  This is done with the private method `recur_step()`.
///
/// 3. Merge the pieces of size `s_i` in order to solve the original problem.  This is done using the private method `merge_step()`.
/// \f[
///   T_{0,N} = T_{0, r} \circ T_{r, s_1} \circ ... \circ T_{..., s_m}.
/// \f]
class GSchurN {
private:
  int N_;  ///< Size of input vector
  int bmod_;  ///< Integer giving binary modulus. Binary pieces of size smaller than `bmod_` are computed using PSchur algorithm.
  double* alpha_; ///< First input polynomial as a function of `acf`.
  double* beta_;  ///< Second input polynomial as a function of `acf`.
  std::vector<int> sbin_;  ///< Vector that records the binary decomposition of `N`.
  int nbin_; ///< Length of `sbin_`.
  GSchur2K** gsb_;  ///< Memory allocation for the binary pieces. There are `sbin_size()` pieces, each of size `sbin_[i]`.
  GSchur2K** gsm_;  ///< Memory allocation for merging the pieces.
  /// Convert integer to modulo-binary representation.
  std::vector<int> int2bin(int n, int b);
  /// Compute `alphan` and `betan` from `alpha0`, `beta0`, `eta0`, and `xi0`.
  void compute_nn(GSchur2K* gsr, int n1, int n2);
  /// Compute `eta2n` and `xi2n` from `etan` and `xin`.
  void compute_2n(GSchur2K* gsr, int n);
  /// Progressive Schur algorithm.
  void prog_step(const double* alpha0, const double* beta0, int n);
  /// Generalized Schur algorithm.
  void recur_step(const double* alpha0, const double* beta0, int si, int layer);
  /// Merge pieces of Schur algorithm contained in `gsb_`.
  void merge_step();
public:
  // double* delta;  ///< The first column of `Toeplitz(acf)^{-1}`.
  // double ldV;   ///< The log-determinant `log(|Toeplitz(acf)|)`.
  /// Constructor.
  GSchurN(int N, int bmod);
  /// Destructor.
  ~GSchurN();
  /// Perform the Generalized Schur algorithm on the input data.
  void compute(double* delta, double& ldV, const double* acf);
};

/// @param[in] N Size of Toeplitz matrix (input to GSchur algorithm is `N-1`).
/// @param[in] bmod Binary modulus (integer).
inline GSchurN::GSchurN(int N, int bmod = 64) {
  N_ = N;
  bmod_ = bmod;
  alpha_ = new double[N_ - 1];
  beta_ = new double[N_ - 1];
  // delta = new double[N_];
  /// sbin_ = (s_1, s_2, ..., s_k)
  sbin_ = int2bin(N_ - 1, bmod_);
  nbin_ = sbin_.size();
  int gs_size = bmod_;
  int gs_layer = log2(ceil((double)sbin_[nbin_ - 1] / bmod_));
  gsb_ = new GSchur2K * [gs_layer + 1];
  gsb_[0] = new GSchur2K(2 * gs_size);
  for (int ii = 0; ii < gs_layer; ++ii) {
    gsb_[ii + 1] = new GSchur2K(2 * gs_size);
    gs_size <<= 1;
  }
  // Length of gsm_ is max(nbin_, 2) - 1
  gsm_ = new GSchur2K * [nbin_ > 2 ? nbin_ : 2 - 1];
  if (nbin_ == 1) {
    gsm_[0] = new GSchur2K(sbin_[0]);
  }
  else {
    int gsm_size = sbin_[0];
    for (int ii = 0; ii < nbin_ - 1; ++ii) {
      gsm_size = gsm_size + sbin_[ii + 1];
      gsm_[ii] = new GSchur2K(gsm_size);
    }
  }
}

inline GSchurN::~GSchurN() {
  delete[] alpha_;
  delete[] beta_;
  // delete[] delta;
  for (int ii = 0; ii <= log2(ceil((double)sbin_[nbin_ - 1] / bmod_)); ++ii) {
    delete gsb_[ii];
  }
  delete[] gsb_;
  if (nbin_ == 1) {
    delete gsm_[0];
  }
  else {
    for (int ii = 0; ii < nbin_ - 1; ++ii) {
      delete gsm_[ii];
    }
  }
  delete[] gsm_;
}

/// Given an integer `n` and modulus `b`, returns the unique integer vector `s = (s_1, ..., s_k)` or `s = (s_0, ..., s_k)` such that:
///
/// - `0 < s_0 < s_1 < ... < s_k`,
/// - `s_i = b * 2^(m_i)` for some integer `m_i`, for `1 <= i <= k`.
/// - If `n` is a multiple of `b`, `s = (s_1, ..., s_k)`.  Otherwise, `s = (s_0, ..., s_k)` and `0 < s_0 < b`.
/// - `n = sum(s)`.
///
/// @param[in] n Integer for whic the modulo-binary represention is computed.
/// @param[in] b Integer giving the binary modulus.
/// @return Integer vector containing the modulo-binary representation.
inline std::vector<int> GSchurN::int2bin(int n, int b = 1) {
  std::vector<int> s;
  int n1 = n / b;
  int dif = n - n1 * b;
  int m = b;
  do {
    if (n1 & 1) {
      s.push_back(m);
    }
    m <<= 1;
  } while ((n1 >>= 1) > 0);
  if (dif) {
    s.insert(s.begin(), dif);
  }
  return s;
}

/// Given polynomials `alpha0(x)`, `beta0(x)`, `eta0(x)`, and `xi0(x)` each represented by a double array `alpha0`, `beta0`, etc., compute polynomials `alphan(x)` and `betan(x)` defined as
///
/// ```
/// alphan(x) = (alpha0(x) * eta0(x) - beta0(x) * xi0(x)) / x^n,
/// betan(x) = (beta0(x) * tilde(eta0(x)) - alpha0(x) * tilde(xi0(x))) / x^n,
/// ```
/// where `tilde(poly(x)) = x^deg(poly) * poly(1/x)`.
///
/// @note By negating the odd parts of the complex conjugate of `fft(eta0, 0, 0, ..., 0)`, we can directly get `fft(0, rev(eta0), 0, ..., 0)`. However, this only works when the amount of zero-padding is also `n`.
///
/// @param[in,out] gsr Pointer to the GSchur storage object containing `alpha0`, `beta0`, etc.  Output is stored in `gsr->alphan` and `grs->betan`.
/// @param[in] n Size of each polynomial `alpha0` and `beta0`, etc.
/// @param[in] npad Integer amount of zero-padding.  The size of each FFT/iFFT is `n + npad`.
inline void GSchurN::compute_nn(GSchur2K* gsr, int n, int npad) {
  int nn = (n + npad) / 2 + 1; // size of elementwise complex products
  // fft of alpha0, beta0, eta0, xi0, etat, xit
  gsr->FFT->fft(gsr->alpha0_fft, gsr->alpha0);
  gsr->FFT->fft(gsr->beta0_fft, gsr->beta0);
  gsr->FFT->fft(gsr->eta0_fft, gsr->eta0);
  gsr->FFT->fft(gsr->xi0_fft, gsr->xi0);
  if (n - npad) {
    for (int ii = 0; ii < n; ++ii) {
      gsr->etat[n - ii] = gsr->eta0[ii];
      gsr->xit[n - ii] = gsr->xi0[ii];
    }
    gsr->FFT->fft(gsr->etat_fft, gsr->etat);
    gsr->FFT->fft(gsr->xit_fft, gsr->xit);
  }
  else {
    for (int ii = 0; ii < n; ++ii) {
      gsr->etat_fft[2 * ii] = std::conj(gsr->eta0_fft[2 * ii]);
      gsr->etat_fft[2 * ii + 1] = -std::conj(gsr->eta0_fft[2 * ii + 1]);
      gsr->xit_fft[2 * ii] = std::conj(gsr->xi0_fft[2 * ii]);
      gsr->xit_fft[2 * ii + 1] = -std::conj(gsr->xi0_fft[2 * ii + 1]);
    }
  }
  // elementwise complex products
  for(int ii=0; ii<nn; ii++) {
    gsr->alphan_fft[ii] = gsr->alpha0_fft[ii] * gsr->eta0_fft[ii] - gsr->xi0_fft[ii] * gsr->beta0_fft[ii];
    gsr->betan_fft[ii] = gsr->beta0_fft[ii] * gsr->etat_fft[ii] - gsr->xit_fft[ii] * gsr->alpha0_fft[ii];
  }
  // complex_mult(gsr->alphan_fft, gsr->alpha0_fft, gsr->eta0_fft, nn);
  // complex_mult_minus(gsr->alphan_fft, gsr->xi0_fft, gsr->beta0_fft, nn);
  // complex_mult(gsr->betan_fft, gsr->beta0_fft, gsr->etat_fft, nn);
  // complex_mult_minus(gsr->betan_fft, gsr->xit_fft, gsr->alpha0_fft, nn);
  // calculate alphan and betan
  gsr->FFT->ifft(gsr->alphan, gsr->alphan_fft);
  gsr->FFT->ifft(gsr->betan, gsr->betan_fft);
}

/// Given polynomials `etan(x)`, `xin(x)`, `eta0(x)`, and `xi0(x)` each represented by a double array `etan`, `xin`, etc., compute polynomials `eta2n(x)` and `xi2n(x)` defined as
///
/// ```
/// eta2n(x) = tilde(xi0(x)) * xin(x) + eta0(x) * etan(x),
/// xi2n(x) = tilde(eta0(x)) * xin(x) + xi0(x) * etan(x).
/// ```
///
/// @param[in,out] gsr Pointer to the GSchur storage object containing `etan`, `xin`, etc.  Output is stored in `gsr->eta2n` and `grs->xin`.
/// @param[in] n Size of each polynomial `etan`, `xin`, etc.
inline void GSchurN::compute_2n(GSchur2K* gsr, int n) {
  int nn = n / 2 + 1; // size of elementwise complex products
  // fft of xin and etan
  gsr->FFT->fft(gsr->xin_fft, gsr->xin);
  gsr->FFT->fft(gsr->etan_fft, gsr->etan);
  // elementwise complex products
  for(int ii=0; ii<nn; ii++) {
    gsr->xi2n_fft[ii] = gsr->etat_fft[ii] * gsr->xin_fft[ii] + gsr->xi0_fft[ii] * gsr->etan_fft[ii];
    gsr->eta2n_fft[ii] = gsr->xit_fft[ii] * gsr->xin_fft[ii] + gsr->eta0_fft[ii] * gsr->etan_fft[ii];
  }
  // complex_mult(gsr->xi2n_fft, gsr->etat_fft, gsr->xin_fft, nn);
  // complex_mult_plus(gsr->xi2n_fft, gsr->xi0_fft, gsr->etan_fft, nn);
  // complex_mult(gsr->eta2n_fft, gsr->xit_fft, gsr->xin_fft, nn);
  // complex_mult_plus(gsr->eta2n_fft, gsr->eta0_fft, gsr->etan_fft, nn);
  // calculate xi2n and eta2n
  gsr->FFT->ifft(gsr->xi2n, gsr->xi2n_fft);
  gsr->FFT->ifft(gsr->eta2n, gsr->eta2n_fft);
}

/// The Progressive Schur (PSchur) algorithm is used for pieces of size less than or equal to the "binary modulus" `bmod_`.  This is because it is faster than the recursive Generalized Schur (GSchur) algorithm for small sizes.  The output of either algorithm is the same, so both can be merged as is done in `recur_step()` and `merge_step()`.  Because of the construction of the GSchur algorithm, PSchur always uses memory in `gsb_[0]`.  The final output of the algorithm is stored in `gsb_[0]->eta2n`, `gsb_[0]->xi2n`, and `gsb_[0]->gamma`.
///
/// @param[in] alpha0 First input real vector.
/// @param[in] beta0 Second input real vector.
/// @param[in] n Size of each input (integer).  This is either `sbin_[0]` or `bmod_`.
inline void GSchurN::prog_step(const double* alpha0, const double* beta0, int n) {
  double *xi1, *xi2, *eta1, *eta2, *swap_ptr; // some pointers to simplify notation
  double alpha1, alpha2;
  // Initialize the memory
  std::fill(gsb_[0]->xin, gsb_[0]->xin + 2 * n, 0);
  std::fill(gsb_[0]->etan, gsb_[0]->etan + 2 * n, 0);
  // assign pointers to members of gsb_[0]
  eta1 = gsb_[0]->etan;
  eta2 = gsb_[0]->etan + n;
  xi1 = gsb_[0]->xin;
  xi2 = gsb_[0]->xin + n;
  // initialize algorithm
  eta1[0] = 1.0;
  xi1[0] = alpha0[0] / beta0[0];
  gsb_[0]->gamma[0] = xi1[0];
  gsb_[0]->beta0[0] = beta0[0] * (1 - xi1[0] * xi1[0]);
  for (int kk = 1; kk < n; ++kk) {
    alpha2 = alpha0[kk];
    gsb_[0]->beta0[kk] = beta0[kk];
    for (int jj = 1; jj <= kk; ++jj) {
      alpha1 = alpha2 - gsb_[0]->gamma[jj - 1] * gsb_[0]->beta0[kk - jj + 1];
      gsb_[0]->beta0[kk - jj + 1] -= gsb_[0]->gamma[jj - 1] * alpha2;
      alpha2 = alpha1;
    }
    gsb_[0]->gamma[kk] = alpha1 / gsb_[0]->beta0[0];
    gsb_[0]->beta0[0] *= 1 - gsb_[0]->gamma[kk] * gsb_[0]->gamma[kk];
    eta2[0] = 1.0;
    xi2[0] = alpha0[0] / beta0[0];
    for (int jj = 1; jj <= kk; ++jj) {
      xi2[jj] = xi1[jj] + gsb_[0]->gamma[kk] * eta1[kk - jj];
      eta2[jj] = eta1[jj] + gsb_[0]->gamma[kk] * xi1[kk - jj];
    }
    swap_ptr = xi1;
    xi1 = xi2;
    xi2 = swap_ptr;
    swap_ptr = eta1;
    eta1 = eta2;
    eta2 = swap_ptr;
  }	
  std::copy(xi1, xi1 + n, gsb_[0]->xi2n);
  std::copy(eta1, eta1 + n, gsb_[0]->eta2n);
  return;
}

/// The Generalized Schur (GSchur) algorithm is a recursive doubling procedure.  It has the same output as PSchur, and shares the same first three inputs as `prog_step()`.  However, it stores the output in `gsb_[layer]` by recursively building up the pieces from `gsb_[0]`.  The final output of the algorithm is stored in `gsb_[layer]->eta2n`, `gsb_[layer]->xi2n`, and `gsb_[layer]->gamma`.
///
/// @param[in] alpha0 First input real vector.
/// @param[in] beta0 Second input real vector.
/// @param[in] n Size of each input (integer).
/// @param[in] layer In which element of `gsb_` array to store computation (integer).
inline void GSchurN::recur_step(const double* alpha0, const double* beta0, int n, int layer) {
  if (n <= bmod_) {
    // for small n use prog_step only
    prog_step(alpha0, beta0, n);
    return;
  }
  prog_step(alpha0, beta0, bmod_);
  int n_ = bmod_;
  for(int mm = 0; mm < layer; mm++) {
    // Prepare the following inputs for compute_nn:
    // alpha0, beta0, xi0, eta0, and gamma
    std::copy(alpha0, alpha0 + 2 * n_, gsb_[mm + 1]->alpha0);
    std::copy(beta0, beta0 + 2 * n_, gsb_[mm + 1]->beta0);
    std::copy(gsb_[mm]->eta2n, gsb_[mm]->eta2n + n_, gsb_[mm + 1]->eta0);
    std::copy(gsb_[mm]->xi2n, gsb_[mm]->xi2n + n_, gsb_[mm + 1]->xi0);
    std::copy(gsb_[mm]->gamma, gsb_[mm]->gamma + n_, gsb_[mm + 1]->gamma);
    // Compute alphan and betan
    compute_nn(gsb_[mm + 1], n_, n_);
    // Compute etan, xin and gamman using recur_step
    recur_step(gsb_[mm + 1]->alphan + n_, gsb_[mm + 1]->betan + n_, n_, mm);
    // Prepare the inputs for compute_2n, including xin and etan
    std::copy(gsb_[mm]->xi2n, gsb_[mm]->xi2n + n_, gsb_[mm + 1]->xin);
    std::fill(gsb_[mm + 1]->xin + n_, gsb_[mm + 1]->xin + 2 * n_, 0);
    std::copy(gsb_[mm]->eta2n, gsb_[mm]->eta2n + n_, gsb_[mm + 1]->etan);
    std::fill(gsb_[mm + 1]->etan + n_, gsb_[mm + 1]->etan + 2 * n_, 0);
    // Compute eta2n and xi2n
    compute_2n(gsb_[mm + 1], 2 * n_);
    // Append gamman after gamma0 to get gamma2n
    std::copy(gsb_[mm]->gamma, gsb_[mm]->gamma + n_, gsb_[mm + 1]->gamma + n_);
    n_ <<= 1;
  }
  return;
}

/// The final output of the algorithm is stored in `gsm_[0]->eta2n`, `gsm_[0]->xi2n`, and `gsm_[0]->gamma`.
inline void GSchurN::merge_step() {
  // calculate each piece of GSchur to a different layer of gsm_
  int layer = log2(ceil((double)sbin_[0] / bmod_));
  recur_step(alpha_, beta_, sbin_[0], layer);
  if (nbin_ == 1) {
    // When vector `sbin_` has only one element, directly pass the results to output.
    std::copy(gsb_[layer]->eta2n, gsb_[layer]->eta2n + sbin_[0], gsm_[0]->eta2n);
    std::copy(gsb_[layer]->xi2n, gsb_[layer]->xi2n + sbin_[0], gsm_[0]->xi2n);
    std::copy(gsb_[layer]->gamma, gsb_[layer]->gamma + sbin_[0], gsm_[0]->gamma);
    return;
  }
  int n_ = 0;
  for(int mm = 0; mm < nbin_ - 1; mm++) {
    n_ += sbin_[mm];
    if(mm) {
      std::copy(gsm_[mm - 1]->eta2n, gsm_[mm - 1]->eta2n + n_, gsm_[mm]->eta0);
      std::copy(gsm_[mm - 1]->xi2n, gsm_[mm - 1]->xi2n + n_, gsm_[mm]->xi0);
    }
    else {	
      std::copy(gsb_[layer]->eta2n, gsb_[layer]->eta2n + n_, gsm_[mm]->eta0);
      std::copy(gsb_[layer]->xi2n, gsb_[layer]->xi2n + n_, gsm_[mm]->xi0);
      std::copy(gsb_[layer]->gamma, gsb_[layer]->gamma + n_, gsm_[nbin_ - 2]->gamma);
    }
    // First paste alpha_{0,n} and beta_{0,n} into gsm_[0]->alpha and
    // gsm_[0]->beta, length is n_
    std::copy(alpha_, alpha_ + n_ + sbin_[mm + 1], gsm_[mm]->alpha0);
    std::copy(beta_, beta_ + n_ + sbin_[mm + 1], gsm_[mm]->beta0);
    // compute_nn uses gsm_[mm]->alpha0, beta0, eta0 and xi0
    // to compute alphan and betan
    compute_nn(gsm_[mm], n_, sbin_[mm + 1]);
    // Computation of recur_step happens inside gsb_[0,1,...,layer].
    // We just need alphan [sbin_[mm]+(1:sbin_[mm+1])]
    layer = log2(ceil((double)sbin_[mm + 1] / bmod_));
    recur_step(gsm_[mm]->alphan + n_, gsm_[mm]->betan + n_, sbin_[mm + 1], layer);
    std::copy(gsb_[layer]->eta2n, gsb_[layer]->eta2n + sbin_[mm + 1], gsm_[mm]->etan);
    std::fill(gsm_[mm]->etan + sbin_[mm + 1], gsm_[mm]->etan + n_ + sbin_[mm + 1], 0);
    std::copy(gsb_[layer]->xi2n, gsb_[layer]->xi2n + sbin_[mm + 1], gsm_[mm]->xin);
    std::fill(gsm_[mm]->xin + sbin_[mm + 1], gsm_[mm]->xin + n_ + sbin_[mm + 1], 0);
    std::copy(gsb_[layer]->gamma, gsb_[layer]->gamma + sbin_[mm + 1], gsm_[nbin_ - 2]->gamma + n_);
    // uses xi_t_fft, xi_0_fft, xi_n_fft (eta), returns xi_0_ifft
    compute_2n(gsm_[mm], sbin_[mm + 1] + n_);
  }
}

/// First generates the input polynomials `alpha0` and `beta0` from the first column `acf` of a symmetric positive-definite Toeplitz matrix `Tz = Toeplitz(acf)`, then applies the Generalized Schur algorithm to each piece in the binary-module representation, and finally merges these pieces to produce the first column of `Tz^{-1}` and calculate the log-determinant of `Tz`. 
///
/// @param[out] delta The first column of `Tz^{-1}`.
/// @param[out] ldV The log-determinant of `Tz`.
/// @param[in] acf The first column of `Tz`.
inline void GSchurN::compute(double* delta, double& ldV, const double* acf) {
  // generate alpha0 and beta0 from acf.
  for (int ii = 0; ii < N_ - 1; ++ii) {
    alpha_[ii] = -1 * acf[ii + 1];
    beta_[ii] = acf[ii];
  }
  // run gschur algorithm
  merge_step();
  // convert gamma into ldV.
  double sigma2 = log(acf[0]);
  ldV = sigma2;
  int k2 = nbin_ > 2 ? nbin_ : 2;
  for (int ii = 0; ii < N_ - 1; ++ii) {
    if (gsm_[k2 - 2]->gamma[ii] < 1) {
      sigma2 += log(1 - gsm_[k2 - 2]->gamma[ii] * gsm_[k2 - 2]->gamma[ii]);
      ldV += sigma2;
    }
  }
  sigma2 = exp(sigma2);
  // convert eta and xi into delta, the first column of inverse Toeplitz matrix.
  std::copy(gsm_[k2 - 2]->eta2n, gsm_[k2 - 2]->eta2n + N_ - 1, delta);
  delta[N_ - 1] = 0;
  delta[0] /= sigma2;
  for (int ii = 1; ii < N_; ++ii) {
    delta[ii] += gsm_[k2 - 2]->xi2n[ii - 1];
    delta[ii] /= sigma2;
  }
  return;
}

#endif

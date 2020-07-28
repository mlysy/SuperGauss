/// @file NormalCirculant.h

#ifndef NormalCirculant_h
#define NormalCirculant_h 1

#include "Circulant.h"

/// @brief The multivariate normal distribution with circulant variance matrix.
///
/// The `N`-dimensional NormalCirculant distribution `z ~ NormalCirculant(uacf)` is defined as
///
/// ```
/// z = (z_1, ..., z_N) ~ Normal(0, Ct = Circulant(uacf)),
/// ```
///
/// where `Circulant(uacf)` denotes a symmetric positive-definite Circulant variance matrix with unique autocorrelation elements `uacf = (uacf_1, ..., uacf_Nu)`, such that the full autocorrelation is
///
/// ```
/// acf = (uacf_1, ..., uacf_Nu, uacf_{Nu-1}, ..., uacf_2),  N even,
///     = (uacf_1, ..., uacf_Nu,   uanc_Nu,   ..., uacf_2),  N odd,
/// ```
///
/// and such that
///
/// ```
/// Ct[i,j] = acf[|i-j|+1],  1 <= i,j <= N.
/// ```
///
/// @todo
///
/// - All FFT outputs contain `Nu` unique (real or complex) elements, so elementwise operations should run through `Nu` and not `N`.  Also means that FFT output sizes can be `Nu` instead of `N`.  Should go fix this in `Circulant.h` and `Toeplitz.h` as well.
class NormalCirculant {
private:
  typedef std::complex<double> dcomplex;
  int N_; ///< Size of multivariate normal.
  int Nu_; ///< Number of unique autocorrelation elements.
  bool Neven_; ///< Whether N is even.
  Circulant* Ct_; ///< Circulant variance matrix.
  double* z_; ///< Circulant density argument.
  double* zsol_; ///< Solution of `Ct_^{-1} z`.
  bool has_z_; ///< Whether z has been set.
  bool has_zsol_; ///< Whether zsol_ has been computed.
  // storage for temporary vectors
  double* vec1_;
  double* vec2_;
  dcomplex* vec1_fft_;
  dcomplex* vec2_fft_;
  RealFFT* rfft_;
  EvenFFT* efft_;
  /// Dot-product between vectors.
  double dot_prod(const double* v1, const double* v2);
public:
  /// Constructor.
  NormalCirculant(int N);
  /// Destructor.
  ~NormalCirculant();
  /// Size of NormalCirculant random vector.
  int size();
  /// Get the NormalCirculant acf.
  void set_acf(const double* uacf);
  /// Get the NormalCirculant acf.
  void get_acf(double* acf);
  /// Check whether the acf has been set.
  bool has_acf();
  /// Set the NormalCirculant z.
  void set_z(const double* z);
  /// Get the NormalCirculant z.
  void get_z(double* z);  
  /// Check whether z has been set.
  bool has_z();
  /// Log-density of NormalCirculant distribution.
  double logdens(const double* z, const double* uacf);
  /// Log-density with internal uacf and z.
  double logdens();
  /// Full gradient of NormalCirculant log-density.
  double grad_full(double* dldz, double* dldu,
		   const double* z, const double* uacf,
		   bool calc_dldz, bool calc_dldu);
  /// Full gradient with internal uacf and z.
  double grad_full(double* dldz, double* dldu,
		   bool calc_dldz, bool calc_dldu);
};

/// @param N Size of NormalCirculant random vector.
inline NormalCirculant::NormalCirculant(int N) {
  N_ = N;
  Nu_ = N/2 + 1;
  Neven_ = (N_%2 == 0);
  Ct_ = new Circulant(N_);
  z_ = new double[N_];
  zsol_ = new double[N_];
  vec1_ = new double[N_];
  vec2_ = new double[N_];
  vec1_fft_ = new dcomplex[N_];
  vec2_fft_ = new dcomplex[N_];
  rfft_ = new RealFFT(N_);
  efft_ = new EvenFFT(N_);
}

inline NormalCirculant::~NormalCirculant() {
  delete Ct_;
  delete[] z_;
  delete[] zsol_;
  delete[] vec1_;
  delete[] vec2_;
  delete[] vec1_fft_;
  delete[] vec2_fft_;
  delete rfft_;
  delete efft_;
}

inline int NormalCirculant::size() {
  return N_;
}

/// @param[in] x First vector of length `N`.
/// @param[in] y Second vector of length `N`.
/// @return Dot product `sum_i=1^n v1_i * v2_i`.
inline double NormalCirculant::dot_prod(const double* v1, const double* v2) {
  double ans = 0.0;
  for (int ii=0; ii<N_; ++ii) {
    ans += v1[ii] * v2[ii];
  }
  return ans;
}

/// @param[in] uacf Unique-elements autocorrelation vector of length `Nu = floor(N/2)+1`.
inline void NormalCirculant::set_acf(const double* uacf) {
  Ct_->set_acf(uacf);
  has_zsol_ = false;
  return;
}

/// @param[out] acf Autocorrelation vector of length `N`.
inline void NormalCirculant::get_acf(double* acf) {
  Ct_->get_acf(acf);
  return;
}

inline bool NormalCirculant::has_acf() { 
  return Ct_->has_acf(); 
}

/// @param[in] z Density argument vector of length `N`.
inline void NormalCirculant::set_z(const double* z) {
  std::copy(z, z + N_, z_);
  has_z_ = true;
  has_zsol_ = false;
  return;
}

/// @param[out] z Density argument vector of length `N`.
inline void NormalCirculant::get_z(double* z) {
  std::copy(z_, z_ + N_, z);
  return;
}

inline bool NormalCirculant::has_z() {
  return has_z_;
}


/// @param[in] z Observation vector of length `N`.
/// @param[in] uacf Unique autocorrelation elements: a vector of length `Nu = floor(N/2)+1`.
/// @return Scalar value of the log-density.
inline double NormalCirculant::logdens(const double* z, const double* uacf) {
  // const double LOG_2PI = 1.837877066409345483560659472811; // log(2pi)
  // double ldens = 0.0;
  // double* Tiz = vec1_;
  // Ct_->set_acf(uacf); // Tz = Toeplitz(acf)
  // Ct_->solve(Tiz, z); // Tiz = Tz^{-1} * z
  // ldens = dot_prod(z, Tiz); // ldens = t(z) * Tz^{-1} * z
  // ldens += Ct_->log_det() + N_ * LOG_2PI;
  // ldens *= -0.5;
  // return ldens;
  Ct_->set_acf(uacf);
  set_z(z);
  return logdens();
}

/// @warning This version will crash if `set_acf()` and `set_z`()` have not been called yet.
inline double NormalCirculant::logdens() {
  const double LOG_2PI = 1.837877066409345483560659472811; // log(2pi)
  double ldens = 0.0;
  // double* Tiz = vec1_;
  // Ct_->set_acf(uacf); // Tz = Toeplitz(acf)
  if(!has_zsol_) Ct_->solve(zsol_, z_); // zsol = Tz^{-1} * z
  ldens = dot_prod(z_, zsol_); // ldens = t(z) * Tz^{-1} * z
  ldens += Ct_->log_det() + N_ * LOG_2PI;
  ldens *= -0.5;
  return ldens;
}

/// Calculates the gradient with respect to each element of `z` and `uacf` of the log-density corresponding to `z ~ NormalCirculant(uacf)`.
///
/// @param[out] dldz Gradient with respect to `z`.  A vector of length `N`.
/// @param[out] dldu Gradient with respect to `uacf`.  A vector of length `Nu = floor(N/2)+1`.
/// @param[in] z Observation vector of length `N`.
/// @param[in] uacf Unique-elements autocorrelation vector of length `Nu`.
/// @param[in] calc_dldz Whether or not to calculate the gradient with respect to `z`.  If `false`, the input vector `dldz` is left unchanged.
/// @param[in] calc_dldu Whether or not to calculate the gradient with respect to `uacf`.  If `false`, the input vector `dldu` is left unchanged.
/// @return The log-density evaluated at `z` and `uacf`.
inline double NormalCirculant::grad_full(double* dldz, double* dldu,
					 const double* z, const double* uacf,
					 bool calc_dldz = true,
					 bool calc_dldu = true) {
  Ct_->set_acf(uacf);
  set_z(z);
  return grad_full(dldz, dldu, calc_dldz, calc_dldu);
  // if(calc_dldz || calc_dldu) {
  //   Ct_->set_acf(uacf);
  //   Ct_->solve(vec1_, z); // vec1_ = Sigma^{-1} z
  // }
  // if(calc_dldz) {
  //   // gradient with respect to z
  //   for (int ii = 0; ii < N_; ii++) {
  //     dldz[ii] = -vec1_[ii];
  //   }
  // }
  // if(calc_dldu) {
  //   // gradient with respect to uacf
  //   // step 1: calculate rho = ifft(fft(dz) * fft(rev(dz)))
  //   std::reverse_copy(vec1_, vec1_ + N_, vec2_); // vec2_ = rev(vec1_)
  //   rfft_->fft(vec1_fft_, vec1_);
  //   rfft_->fft(vec2_fft_, vec2_);
  //   for(int ii=0; ii<N_; ii++) {
  //     vec1_fft_[ii] *= vec2_fft_[ii];
  //   }
  //   rfft_->ifft(vec1_, vec1_fft_); // vec1_ = rev(rho)
  //   // step 2: calculate kappa = N * ifft(1/psd) = fft(1/psd)
  //   Ct_->get_psd(vec2_);
  //   for(int ii=0; ii<Nu_; ii++) {
  //     vec2_[ii] = 1.0/vec2_[ii];
  //   }
  //   efft_->fft(vec2_, vec2_);
  //   // step 3: combine results
  //   for(int ii=0; ii<Nu_; ii++) {
  //     dldu[ii] = vec1_[N_-ii-1] - vec2_[ii];
  //   }
  //   dldu[0] *= 0.5;
  //   if(Neven_) dldu[Nu_-1] *= 0.5;
  // }
  // return;
}


/// @warning This version will crash if `set_acf()` and `set_z()` have not been called yet.
inline double NormalCirculant::grad_full(double* dldz, double* dldu,
					 bool calc_dldz = true,
					 bool calc_dldu = true) {
  // if(calc_dldz || calc_dldu) {
  //   Ct_->set_acf(uacf);
  //   Ct_->solve(vec1_, z); // vec1_ = Sigma^{-1} z
  // }
  if(!has_zsol_) Ct_->solve(zsol_, z_);
  if(calc_dldz) {
    // gradient with respect to z
    for (int ii = 0; ii < N_; ii++) {
      dldz[ii] = -zsol_[ii];
    }
  }
  if(calc_dldu) {
    // gradient with respect to uacf
    // step 1: calculate rho = ifft(fft(dz) * fft(rev(dz)))
    std::reverse_copy(zsol_, zsol_ + N_, vec2_); // vec2_ = rev(vec1_)
    rfft_->fft(vec1_fft_, zsol_);
    rfft_->fft(vec2_fft_, vec2_);
    for(int ii=0; ii<N_; ii++) {
      vec1_fft_[ii] *= vec2_fft_[ii];
    }
    rfft_->ifft(vec1_, vec1_fft_); // vec1_ = rev(rho)
    // step 2: calculate kappa = N * ifft(1/psd) = fft(1/psd)
    Ct_->get_psd(vec2_);
    for(int ii=0; ii<Nu_; ii++) {
      vec2_[ii] = 1.0/vec2_[ii];
    }
    efft_->fft(vec2_, vec2_);
    // step 3: combine results
    for(int ii=0; ii<Nu_; ii++) {
      dldu[ii] = vec1_[N_-ii-1] - vec2_[ii];
    }
    dldu[0] *= 0.5;
    if(Neven_) dldu[Nu_-1] *= 0.5;
  }
  return logdens();
}


#endif

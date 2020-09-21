/// @file NormalToeplitz.h

#ifndef NormalToeplitz_h
#define NormalToeplitz_h 1

#include "Toeplitz.h"

/// @brief The multivariate normal distribution with Toeplitz variance matrix.
///
/// The NormalToeplitz distribution `z ~ NormalToeplitz(acf)` is defined as
///
/// ```
/// z = (z_1, ..., z_N) ~ Normal(0, Tz = Toeplitz(acf)),
/// ```
///
/// where `Toeplitz(acf)` denotes a symmetric positive-definite Toeplitz variance matrix with autocorrelation (i.e., first row/column) `acf = (acf_1, ..., acf_N)`, i.e.,
///
/// ```
/// Tz[i,j] = acf[|i-j|+1],  1 <= i,j <= N.
/// ```
///
class NormalToeplitz {
private:
  int N_; ///< Size of multivariate normal.
  Toeplitz* Tz_; ///< Toeplitz variance matrix.
  double* z_; ///< Toeplitz density argument.
  double* zsol_; ///< Solution of `Tz_^{-1} z`.
  bool has_z_; ///< Whether z has been set.
  bool has_zsol_; ///< Whether zsol_ has been computed.
  // storages for temporary vectors
  double* vec2;
  double* vec3;
  double* vec4;
  double* phi;
  /// Dot-product between vectors.
  double dot_prod(const double* v1, const double* v2);
public:
  /// Constructor.
  NormalToeplitz(int N);
  /// Destructor.
  ~NormalToeplitz();
  /// Size of NormalToeplitz random vector.
  int size();
  /// Set the NormalToeplitz acf.
  void set_acf(const double* acf);
  /// Get the NormalToeplitz acf.
  void get_acf(double* acf);  
  /// Check whether the acf has been set.
  bool has_acf();
  /// Set the NormalToeplitz z.
  void set_z(const double* z);
  /// Get the NormalToeplitz z.
  void get_z(double* z);  
  /// Check whether z has been set.
  bool has_z();
  /// Log-density of NormalToeplitz distribution.
  double logdens(const double* z, const double* acf);
  /// Log-density with internal acf and z.
  double logdens();
  /// Gradient of NormalToeplitz log-density.
  void grad(double* dldt, const double* z, const double* dzdt,
	    const double* acf, const double* dadt, int n_theta);
  /// Gradient with internal acf and z.
  void grad(double* dldt, const double* dzdt,
	    const double* dadt, int n_theta);
  /// Hessian of NormalToeplitz log-density.
  void hess(double* d2ldt,
	    const double* z, const double* dzdt, const double* d2zdt,
	    const double* acf, const double* dadt, const double* d2adt,
	    int n_theta);
  /// Hessian with internal acf and z.
  void hess(double* d2ldt,
	    const double* dzdt, const double* d2zdt,
	    const double* dadt, const double* d2adt,
	    int n_theta);
  /// Full gradient of NormalToeplitz log-density.
  double grad_full(double* dldz, double* dlda,
		   const double* z, const double* acf,
		   bool calc_dldz, bool calc_dlda);
  /// Full gradient with internal acf and z.
  double grad_full(double* dldz, double* dlda,
		   bool calc_dldz, bool calc_dlda);
};

/// @param N Size of NormalToeplitz random vector.
inline NormalToeplitz::NormalToeplitz(int N) {
  N_ = N;
  Tz_ = new Toeplitz(N_);
  z_ = new double[N_];
  zsol_ = new double[N_];
  vec2 = new double[N_];
  vec3 = new double[N_];
  vec4 = new double[N_];
  phi = new double[N_];	
}

inline NormalToeplitz::~NormalToeplitz() {
  delete Tz_;
  delete[] z_;
  delete[] zsol_;
  delete[] vec2;
  delete[] vec3;
  delete[] vec4;
  delete[] phi;
}

inline int NormalToeplitz::size() {
  return N_;
}

/// @param[in] v1 First vector of length `N`.
/// @param[in] v2 Second vector of length `N`.
/// @return Dot product `sum_i=1^n v1_i * v2_i`.
inline double NormalToeplitz::dot_prod(const double* v1, const double* v2) {
  double ans = 0.0;
  for (int ii=0; ii<N_; ++ii) {
    ans += v1[ii] * v2[ii];
  }
  return ans;
}

/// @param[in] acf Autocorrelation vector of length `N`.
///
/// @note For `logdens()`, `grad()`, `grad_full`, and `hess()` methods with multiple `z` but common `acf`, the expensive call to `Toeplitz::solve()` needn't be recomputed.
inline void NormalToeplitz::set_acf(const double* acf) {
  Tz_->set_acf(acf);
  has_zsol_ = false;
  return;
}

/// @param[in] acf Autocorrelation vector of length `N`.
inline void NormalToeplitz::get_acf(double* acf) {
  Tz_->get_acf(acf);
  return;
}

inline bool NormalToeplitz::has_acf() { 
  return Tz_->has_acf(); 
}

/// @param[in] z Density argument vector of length `N`.
inline void NormalToeplitz::set_z(const double* z) {
  std::copy(z, z + N_, z_);
  has_z_ = true;
  has_zsol_ = false;
  return;
}

/// @param[out] z Density argument vector of length `N`.
inline void NormalToeplitz::get_z(double* z) {
  std::copy(z_, z_ + N_, z);
  return;
}

inline bool NormalToeplitz::has_z() {
  return has_z_;
}

/// @param[in] z Observation vector of length `N`.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @return Scalar value of the log-density.
inline double NormalToeplitz::logdens(const double* z, const double* acf) {
  Tz_->set_acf(acf); // Tz = Toeplitz(acf)
  set_z(z);
  return logdens();
}

/// @warning This version will crash if `set_acf()` and `set_z`()` have not been called yet.
inline double NormalToeplitz::logdens() {
  const double LOG_2PI = 1.837877066409345483560659472811; // log(2pi)
  double ldens = 0.0;
  // Tz_->set_acf(acf); // Tz = Toeplitz(acf)
  if(!has_zsol_) Tz_->solve(zsol_, z_); // zsol = Tz^{-1} * z
  ldens = dot_prod(z_, zsol_); // ldens = t(z) * Tz^{-1} * z
  ldens += Tz_->log_det() + N_ * LOG_2PI;
  ldens *= -0.5;
  return ldens;
}

/// Calculates the gradient with respect to `theta` of the loglikelihood corresponding to
///
/// ```
/// z(theta) ~ NormalToeplitz(acf(theta)).
/// ```
///
/// @param[out] dldt Gradient of the loglikelihood.  A vector of length `n_theta`.
/// @param[in] z Observation vector of length `N`.
/// @param[in] dzdt Gradient of `z` with respect to `theta`.  A vector of length `N * n_theta` corresponding to the gradient matrix of size `N x n_theta` flattened in column-major order.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[in] dadt Gradient of `acf` with respect to `theta`.  A vector of length `N * n_theta` corresponding to the gradient matrix of size `N x n_theta` flattened in column-major order.
inline void NormalToeplitz::grad(double* dldt,
				 const double* z, const double* dzdt, 
				 const double* acf, const double* dadt,
				 int n_theta) {
  Tz_->set_acf(acf); // Tz = Toeplitz(acf)
  set_z(z);
  grad(dldt, dzdt, dadt, n_theta);
}

/// @warning In this version `dzdt` and `dadt` must still be supplied.  This version will crash if `set_acf()` and `set_z`()` have not been called yet.
inline void NormalToeplitz::grad(double* dldt,
				 const double* dzdt, 
				 const double* dadt,
				 int n_theta) {
  // Tz_->set_acf(acf); // Tz = Toeplitz(acf)
  if(!has_zsol_) Tz_->solve(zsol_, z_); // zsol = Tz^{-1} * z
  for(int ii = 0; ii < n_theta; ++ii) {
    Tz_->prod(vec2, zsol_, &dadt[ii * N_]);
    dldt[ii] = .5 * dot_prod(vec2, zsol_);
    dldt[ii] -= dot_prod(&dzdt[ii * N_], zsol_);
    dldt[ii] -= .5 * Tz_->trace_grad(&dadt[ii * N_]);
  }
  return;
}



/// Calculates the Hessian matrix with respect to `theta` of the loglikelihood corresponding to
///
/// ```
/// z(theta) ~ NormalToeplitz(acf(theta)).
/// ```
///
/// @param[out] d2ldt Hessian of the loglikelihood.  A vector of length `n_theta * n_theta` corresponding to the Hessian matrix of size `n_theta x n_theta` flattened in column-major order.
/// @param[in] z Observation vector of length `N`.
/// @param[in] dzdt Gradient of `z` with respect to `theta`.  A vector of length `N * n_theta` corresponding to the gradient matrix of size `N x n_theta` flattened in column-major order.
/// @param[in] d2zdt Hessian of `z` with respect to `theta`.  A vector of length `N * n_theta * n_theta` corresponding to the Hessian tensor of size `N x n_theta x n_theta` flattened in column-major order (i.e., leftmost dimension running fastest).
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[in] dadt Gradient of `acf` with respect to `theta`.  A vector of length `N * n_theta` corresponding to the gradient matrix of size `N x n_theta` flattened in column-major order.
/// @param[in] d2adt Hessian of `acf` with respect to `theta`.  A vector of length `N * n_theta * n_theta` corresponding to the Hessian tensor of size `N x n_theta x n_theta` flattened in column-major order.
inline void NormalToeplitz::hess(double* d2ldt,
				 const double* z,
				 const double* dzdt,
				 const double* d2zdt,
				 const double* acf,
				 const double* dadt,
				 const double* d2adt,
				 int n_theta) {
  Tz_->set_acf(acf);
  set_z(z);
  hess(d2ldt, dzdt, d2zdt, dadt, d2adt, n_theta);
  return;
}

/// @warning In this version `dzdt`, `d2zdt`, `dadt`, and `d2adt` must still be supplied.  This version will crash if `set_acf()` and `set_z()` have not been called yet.
inline void NormalToeplitz::hess(double* d2ldt,
				 const double* dzdt,
				 const double* d2zdt,
				 const double* dadt,
				 const double* d2adt,
				 int n_theta) {
  // Tz_->set_acf(acf);
  if(!has_zsol_) Tz_->solve(zsol_, z_); // zsol = Tz^{-1} * z
  double ans;
  std::fill(d2ldt, d2ldt + n_theta * n_theta, 0.0);
  for(int ii = 0; ii < n_theta; ++ii) {
    for(int jj = 0; jj <= ii; ++jj) {
      Tz_->prod(vec4, zsol_, &dadt[jj * N_]);
      Tz_->prod(vec3, zsol_, &dadt[ii * N_]);
      ans = dot_prod(&d2zdt[(ii * n_theta + jj) * N_], zsol_);
      Tz_->solve(vec2, vec4);
      ans -= dot_prod(&dzdt[ii * N_], vec2);
      ans += dot_prod(vec3, vec2);
      Tz_->solve(vec2, vec3);
      ans -= dot_prod(&dzdt[jj * N_], vec2);
      Tz_->solve(vec2, &dzdt[jj * N_]);
      ans += dot_prod(&dzdt[ii * N_], vec2);
      ans *= 2.0;
      Tz_->prod(vec2, zsol_, &d2adt[(ii * n_theta + jj) * N_]);
      ans -= dot_prod(zsol_, vec2);
      ans += Tz_->trace_grad(&d2adt[(ii * n_theta + jj) * N_]);
      ans -= Tz_->trace_hess(&dadt[ii * N_], &dadt[jj * N_]);
      d2ldt[ii * n_theta + jj] = -.5 * ans;
    }
  }
  if(n_theta > 1) {
    // copy other triangular half of hessian
    for(int ii = 0; ii < n_theta; ++ii) {
      for(int jj = ii + 1; jj < n_theta; ++jj) {
	d2ldt[ii * n_theta + jj] = d2ldt[jj * n_theta + ii];
      }
    }
  }
  return;
}


/// Calculates the gradient with respect to each element of `z` and `acf` of the log-density corresponding to `z ~ NormalToeplitz(acf)`.
///
/// @param[out] dldz Gradient with respect to `z`.  A vector of length `N`.
/// @param[out] dlda Gradient with respect to `acf`.  A vector of length `N`.
/// @param[in] z Observation vector of length `N`.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[in] calc_dldz Whether or not to calculate the gradient with respect to `z`.  If `false`, the input vector `dldz` is left unchanged.
/// @param[in] calc_dlda Whether or not to calculate the gradient with respect to `acf`.  If `false`, the input vector `dlda` is left unchanged.
/// @return The log-density evaluated at `z` and `acf`.
inline double NormalToeplitz::grad_full(double* dldz, double* dlda,
					const double* z, const double* acf,
					bool calc_dldz = true,
					bool calc_dlda = true) {
  // if(calc_dldz || calc_dlda) {
  //   Tz_->set_acf(acf);
  // }
  Tz_->set_acf(acf);
  set_z(z);
  return grad_full(dldz, dlda, calc_dldz, calc_dlda);
}

/// @warning This version will crash if `set_acf()` and `set_z()` have not been called yet.
inline double NormalToeplitz::grad_full(double* dldz, double* dlda,
					bool calc_dldz = true,
					bool calc_dlda = true) {
  // if(calc_dldz || calc_dlda) {
  //   // Tz_->set_acf(acf);
  //   Tz_->solve(zsol_, z);	
  // }
  if(!has_zsol_) Tz_->solve(zsol_, z_);
  if(calc_dldz) {
    // gradient with respect to z
    for (int ii = 0; ii < N_; ii++) {
      dldz[ii] = -zsol_[ii];
    }
  }
  if(calc_dlda) {
    // gradient with respect to acf
    vec2[0] = 1.0;
    std::fill(vec2 + 1, vec2 + N_, 0.0);
    Tz_->solve(vec3, vec2);
    // dlda = upper.toep(Vz) %*% Vz = ip
    double tau1 = vec3[0];
    std::fill(phi, phi + N_, 0.0);
    phi[0] = zsol_[0];
    Tz_->prod(dlda, zsol_, phi, zsol_);
    // vec2 = (N_:1 * tau)
    vec4[0] = 0.0;
    for (int ii = 1; ii < N_; ++ii) {
      vec4[ii] = vec3[N_ - ii];
    }
    for (int ii = 0; ii < N_; ++ii) {
      vec2[ii] = (N_ - ii) * vec3[ii];
    }
    // vec3 = upper.toep(tau) %*% (N_:1 * tau) = tr
    phi[0] = vec3[0];
    Tz_->prod(vec3, vec2, phi, vec3);
    // vec2 = (N_:1 * tau2)
    for (int ii = 0; ii < N_; ++ii) {
      vec2[ii] = (N_ - ii) * vec4[ii];
    } 
    // vec4 = upper.toep(tau2) %*% (N_:1 * tau2)
    phi[0] = vec4[0];
    Tz_->prod(vec4, vec2, phi, vec4);
    // vec3 = (vec3 - vec4) / tau[1] = tr, dlda = ip - tr
    for (int ii = 0; ii < N_; ++ii) {
      vec3[ii] -= vec4[ii];
      vec3[ii] /= tau1;
      dlda[ii] -= vec3[ii];
    }
    dlda[0] *= .5;
  }
  return logdens();
}

#endif

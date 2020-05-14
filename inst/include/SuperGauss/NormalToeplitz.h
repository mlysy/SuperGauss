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
  // storages for temporary vectors
  double* vec1;
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
  /// Log-density of NormalToeplitz distribution.
  double logdens(const double* z, const double* acf);
  /// Gradient of NormalToeplitz loglikelihood.
  void grad(double* dldt, const double* z, const double* dzdt,
	    const double* acf, const double* dadt, int n_theta);
  /// Hessian of NormalToeplitz loglikelihood.
  void hess(double* d2ldt,
	    const double* z, const double* dzdt, const double* d2zdt,
	    const double* acf, const double* dadt, const double* d2adt,
	    int n_theta);
  /// Full gradient of NormalToeplitz log-density.
  void grad_full(double* dldz, double* dlda,
		 const double* z, const double* acf,
		 bool calc_dldz, bool calc_dlda);
};

/// @param N Size of NormalToeplitz random vector.
inline NormalToeplitz::NormalToeplitz(int N) {
  N_ = N;
  Tz_ = new Toeplitz(N_);
  vec1 = new double[N_];
  vec2 = new double[N_];
  vec3 = new double[N_];
  vec4 = new double[N_];
  phi = new double[N_];	
}

inline NormalToeplitz::~NormalToeplitz() {
  delete Tz_;
  delete[] vec1;
  delete[] vec2;
  delete[] vec3;
  delete[] vec4;
  delete[] phi;
}

inline int NormalToeplitz::size() {
  return N_;
}

/// @param[in] x First vector of length `N`.
/// @param[in] y Second vector of length `N`.
/// @return Dot product `sum_i=1^n v1_i * v2_i`.
inline double NormalToeplitz::dot_prod(const double* v1, const double* v2) {
  double ans = 0.0;
  for (int ii=0; ii<N_; ++ii) {
    ans += v1[ii] * v2[ii];
  }
  return ans;
}


/// @param[in] z Observation vector of length `N`.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[out] Scalar value of the log-density.
inline double NormalToeplitz::logdens(const double* z, const double* acf) {
  const double LOG_2PI = 1.837877066409345483560659472811; // log(2pi)
  double ldens = 0.0;
  Tz_->set_acf(acf); // Tz = Toeplitz(acf)
  Tz_->solve(vec1, z); // vec1 = Tz^{-1} * z
  ldens = dot_prod(z, vec1); // ldens = t(z) * Tz^{-1} * z
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
  Tz_->solve(vec1, z); // vec1 = Tz^{-1} * z
  for(int ii = 0; ii < n_theta; ++ii) {
    Tz_->prod(vec2, vec1, &dadt[ii * N_]);
    dldt[ii] = .5 * dot_prod(vec1, vec2);
    dldt[ii] -= dot_prod(&dzdt[ii * N_], vec1);
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
  Tz_->solve(vec1, z);
  double ans;
  std::fill(d2ldt, d2ldt + n_theta * n_theta, 0.0);
  for(int ii = 0; ii < n_theta; ++ii) {
    for(int jj = 0; jj <= ii; ++jj) {
      Tz_->prod(vec4, vec1, &dadt[jj * N_]);
      Tz_->prod(vec3, vec1, &dadt[ii * N_]);
      ans = dot_prod(&d2zdt[(ii * n_theta + jj) * N_], vec1);
      Tz_->solve(vec2, vec4);
      ans -= dot_prod(&dzdt[ii * N_], vec2);
      ans += dot_prod(vec3, vec2);
      Tz_->solve(vec2, vec3);
      ans -= dot_prod(&dzdt[jj * N_], vec2);
      Tz_->solve(vec2, &dzdt[jj * N_]);
      ans += dot_prod(&dzdt[ii * N_], vec2);
      ans *= 2.0;
      Tz_->prod(vec2, vec1, &d2adt[(ii * n_theta + jj) * N_]);
      ans -= dot_prod(vec1, vec2);
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
inline void NormalToeplitz::grad_full(double* dldz, double* dlda,
				      const double* z, const double* acf,
				      bool calc_dldz = true,
				      bool calc_dlda = true) {
  if(calc_dldz || calc_dlda) {
    Tz_->set_acf(acf);
    Tz_->solve(vec1, z);	
  }
  if(calc_dldz) {
    // gradient with respect to z
    // Tz_->solve(dldz, z);
    for (int ii = 0; ii < N_; ii++) {
      // dldz[ii] = -dldz[ii];
      dldz[ii] = -vec1[ii];
    }
  }
  if(calc_dlda) {
    // gradient with respect to acf
    // Tz_->solve(vec1, z);	
    vec2[0] = 1.0;
    std::fill(vec2 + 1, vec2 + N_, 0.0);
    Tz_->solve(vec3, vec2);
    // dlda = upper.toep(Vz) %*% Vz = ip
    double tau1 = vec3[0];
    std::fill(phi, phi + N_, 0.0);
    phi[0] = vec1[0];
    Tz_->prod(dlda, vec1, phi, vec1);
    // vec2 = (N_:1 * tau)
    vec4[0] = 0.0;
    for (int ii = 1; ii < N_; ++ii) {
      vec4[ii] = vec3[N_ - ii];
    }
    for (int ii = 0; ii < N_; ++ii) {
      vec2[ii] = (N_ - ii) * vec3[ii];
    }
    // vec1 = upper.toep(tau) %*% (N_:1 * tau) = tr
    phi[0] = vec3[0];
    Tz_->prod(vec1, vec2, phi, vec3); 
    // vec2 = (N_:1 * tau2)
    for (int ii = 0; ii < N_; ++ii) {
      vec2[ii] = (N_ - ii) * vec4[ii];
    } 
    // vec3 = upper.toep(tau2) %*% (N_:1 * tau2)
    phi[0] = vec4[0];
    Tz_->prod(vec3, vec2, phi, vec4); 
    // vec1 = (vec1 - vec3) / tau[1] = tr, dlda = ip - tr
    for (int ii = 0; ii < N_; ++ii) {
      vec1[ii] -= vec3[ii];
      vec1[ii] /= tau1;
      dlda[ii] -= vec1[ii];
    }
    dlda[0] *= .5;
  }
  return;
}


/*

//////////////////////////////// overloaded version //////////////////////////////////////////

/// Overloaded Log-Density.
///
/// @param[in] z length-N vector of observation.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
/// @param[out] double number of the log-density.
inline double NormalToeplitz::logdens(double* z, Toeplitz* Tz) {
double ldens = 0;
// acf is already stored in Tz
Tz->solveVec(vec1, z); // vec1 = Tz^{-1} * z
ldens = dot_prod(z, vec1, N_); // ldens = t(z) * Tz^{-1} * z
ldens += Tz->log_det() + N_ * log2Pi;
ldens *= -0.5;
return ldens;
}

/// Overloaded Log-likelihood Gradient.
///
/// @param[out] dldt length-p vector of the gradient.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
/// @param[in] dadt length-N*p vector of observation derivatives, elememts dadt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
inline void NormalToeplitz::grad(double* dldt, double* z, double* dzdt,
Toeplitz* Tz, double* dadt) {
// acf is already stored in Tz
Tz->solveVec(vec1, z); // vec1 = Tz^{-1} * z
double tmp;
for (int ii = 0; ii < n_theta; ++ii) {
T2_->set_acf(&dadt[ii * N_]);
T2_->multVec(vec2, vec1);
dldt[ii] = dot_prod(vec1, vec2, N_) / 2;
dldt[ii] -= dot_prod(&dzdt[ii * N_], vec1, N_);
dldt[ii] -= Tz->traceProd(&dadt[ii * N_]) / 2;
}
}

/// Overloaded function for the Hessian matrix of Log-likelihood.
///
/// @param[out] d2ldt length-p*p vector of the Hessian, elememts d2ldt[0:p + ii*p] is the (ii+1)-th column of Hessian matrix.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] d2zdt length-N*p*p vector of observation derivatives, elememts d2zdt[0:N + (ii*p+jj)*N] is the second derivative of z with respect to the ii-th parameter and then the jj-th parameter.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
/// @param[in] dadt length-N*p vector of observation derivatives, elememts dadt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
/// @param[in] d2adt length-N*p*p vector of observation derivatives, elememts d2adt[0:N + (ii*p+jj)*N] is the second derivative of acf with respect to the ii-th parameter and then the jj-th parameter.
inline void NormalToeplitz::hess(double* d2ldt,
double* z, double* dzdt, double* d2zdt,
Toeplitz* Tz, double* dadt, double* d2adt) {
// acf is already stored in Tz
Tz->solveVec(vec1, z);
double ans;
for (int ii = 0; ii < n_theta; ++ii) {
for (int jj = 0; jj <= ii; ++jj) {
T2_->set_acf(&dadt[jj * N_]);
T2_->multVec(vec4, vec1);
T2_->set_acf(&dadt[ii * N_]);
T2_->multVec(vec3, vec1);

ans = dot_prod(&d2zdt[(ii * n_theta + jj) * N_], vec1, N_);
// temp0 = solve(Tz, temp2)
Tz->solveVec(vec2, vec4);
ans -= dot_prod(&dzdt[ii * N_], vec2, N_);
ans += dot_prod(vec3, vec2, N_);
// temp0 = solve(Tz, temp1)
Tz->solveVec(vec2, vec3);
ans -= dot_prod(&dzdt[jj * N_], vec2, N_);
// temp0 = solve(Tz, dZ[,jj])
Tz->solveVec(vec2, &dzdt[jj * N_]);
ans += dot_prod(&dzdt[ii * N_], vec2, N_);
ans *= 2;

T2_->set_acf(&d2adt[(ii * n_theta + jj) * N_]);
T2_->multVec(vec2, vec1);
ans -= dot_prod(vec1, vec2, N_);
ans += Tz->traceProd(&d2adt[(ii * n_theta + jj) * N_]);
ans -= Tz->traceDeriv(&dadt[ii * N_], &dadt[jj * N_]);

d2ldt[ii * n_theta + jj] = -ans / 2;
}
}

if (n_theta > 1) {
for (int ii = 0; ii < n_theta; ++ii) {
for (int jj = ii + 1; jj < n_theta; ++jj) {
d2ldt[ii * n_theta + jj] = d2ldt[jj * n_theta + ii];
}
}
}

}

/// Overloaded function for the Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
///
/// @param[out] dldz length-N vector of the log-density with respect to the observation z.
/// @param[out] dlda length-N vector of the log-density with respect to the auto-covariance acf.
/// @param[in] z length-N vector of observation.
/// @param[in] Tz size-N Toeplitz object with acf inputed.
inline void NormalToeplitz::grad_full(double* dldz, double* dlda,
double* z, Toeplitz* Tz) {
// acf is already stored in Tz

// gradient with respect to z
Tz->solveVec(dldz, z);
for (int ii = 0; ii < N_; ii++) {
dldz[ii] = -dldz[ii];
}

// gradient with respect to acf
Tz->solveVec(vec1, z); // vec1 = Vz

vec2[0] = 1;
std::fill(vec2 + 1, vec2 + N_, 0); // vec2 = [1,0,0,...,0]
Tz->solveVec(vec3, vec2); // vec3 = tau
double tau1 = vec3[0];

T2_->set_acf(vec1);
T2_->mult0Vec(dlda, vec1); // dlda = upper.toep(Vz) %*% Vz = ip

vec4[0] = 0;
for (int ii = 1; ii < N_; ++ii) {
vec4[ii] = vec3[N_ - ii];
} // vec4 = tau2

for (int ii = 0; ii < N_; ++ii) {
vec2[ii] = (N_ - ii) * vec3[ii];
} // vec2 = (N_:1 * tau)

T2_->set_acf(vec3);
T2_->mult0Vec(vec1, vec2); // vec1 = upper.toep(tau) %*% (N_:1 * tau) = tr

for (int ii = 0; ii < N_; ++ii) {
vec2[ii] = (N_ - ii) * vec4[ii];
} // vec2 = (N_:1 * tau2)
T2_->set_acf(vec4);
T2_->mult0Vec(vec3, vec2); // vec3 = upper.toep(tau2) %*% (N_:1 * tau2)

for (int ii = 0; ii < N_; ++ii) {
vec1[ii] -= vec3[ii];
vec1[ii] /= tau1;
dlda[ii] -= vec1[ii];
} // vec1 = (vec1 - vec3) / tau[1] = tr
// dlda = ip - tr

dlda[0] /= 2;
}

*/

#endif

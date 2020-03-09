/// @file Toeplitz.h

///////////////////////////////////////////////
// Toeplitz matrix class
///////////////////////////////////////////////

#ifndef NormalToeplitz_h
#define NormalToeplitz_h 1

#include "Toeplitz.h"
# define log2Pi 1.83787706641L // log of 2 PI value


/// Class for Computation involving Toeplitz matrix.
///
/// Model: z ~ N(0, V), where V is a Toeplitz matrix with parameter theta
class NormalToeplitz {
private:
  int N_; ///< Size of multivariate normal.
  Toeplitz* Tz_; ///< Toeplitz variance matrix.
  int p_; // number of parameters theta	
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
  NormalToeplitz(int N, int p);
  /// Destructor.
  ~NormalToeplitz();

  /// Size
  int size();

  /// Dimension
  int dim();

  /// Log-Density.
  double logdens(const double* z, const double* acf);

  /// Log-likelihood Gradient.
  void grad(double* dldt, const double* z, const double* dzdt,
	    const double* acf, const double* dadt);

  /// Hessian matrix of Log-likelihood.
  void hess(double* d2ldt,
	    const double* z, const double* dzdt, const double* d2zdt,
	    const double* acf, const double* dadt, const double* d2adt);

  /// Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
  void grad_full(double* dldz, double* dlda,
		 const double* z, const double* acf);

  /*
 /// Overloaded Log-Density.
 double logdens(double* z, Toeplitz* Tz);

 /// Overloaded Log-likelihood Gradient.
 void grad(double* dldt,
 double* z, double* dzdt,
 Toeplitz* Tz, double* dacfdt);

 /// Overloaded function for the Hessian matrix of Log-likelihood.
 void hess(double* d2ldt,
 double* z, double* dzdt, double* d2zdt,
 Toeplitz* Tz, double* dacfdt, double* d2adt);

 /// Overloaded function for the Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
 void grad_full(double* dldz, double* dlda,
 double* z, Toeplitz* Tz);
  */
};

/// Constructor.
///
/// @param N length of observation.
/// @param p_ number of parameters.
/// @param hasToep whether Toeplitz object Tz is imported from outside or generated within the class.
inline NormalToeplitz::NormalToeplitz(int N, int p) {
  N_ = N;
  p_ = p;
  Tz_ = new Toeplitz(N_);
  vec1 = new double[N_];
  vec2 = new double[N_];
  vec3 = new double[N_];
  vec4 = new double[N_];
  phi = new double[N_];	
}

/// Destructor
inline NormalToeplitz::~NormalToeplitz() {
  delete Tz_;
  delete vec1;
  delete vec2;
  delete vec3;
  delete vec4;
  delete phi;
}

inline int NormalToeplitz::size() {
  return N_;
}

inline int NormalToeplitz::dim() {
  return p_;
}

/// @param[in] x First vector of length `N`.
/// @param[in] y Second vector of length `N`.
/// @param[out] Dot product `sum_i=1^n v1_i * v2_i.
inline double NormalToeplitz::dot_prod(const double* v1, const double* v2) {
  double ans = 0.0;
  for (int ii=0; ii<N_; ++ii) {
    ans += v1[ii] * v2[ii];
  }
  return ans;
}


/// Log-Density.
///
/// @param[in] z length-N vector of observation.
/// @param[in] acf length-N vector of auto-covariance.
/// @param[out] double number of the log-density.
inline double NormalToeplitz::logdens(const double* z, const double* acf) {
  double ldens = 0;
  Tz_->setAcf(acf); // Tz = Toeplitz(acf)
  Tz_->solve(vec1, z); // vec1 = Tz^{-1} * z
  ldens = dot_prod(z, vec1); // ldens = t(z) * Tz^{-1} * z
  ldens += Tz_->logDet() + N_ * log2Pi;
  ldens *= -0.5;
  return ldens;
}


/// Log-likelihood Gradient.
///
/// @param[out] dldt length-p vector of the gradient.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] acf length-N vector of auto-covariance.
/// @param[in] dadt length-N*p vector of observation derivatives, elememts dadt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
inline void NormalToeplitz::grad(double* dldt,
				 const double* z, const double* dzdt, 
				 const double* acf, const double* dadt) {
  Tz_->setAcf(acf); // Tz = Toeplitz(acf)
  Tz_->solve(vec1, z); // vec1 = Tz^{-1} * z
  for (int ii = 0; ii < p_; ++ii) {
    Tz_->product(vec2, vec1, &dadt[ii * N_]);
    dldt[ii] = dot_prod(vec1, vec2) / 2;
    dldt[ii] -= dot_prod(&dzdt[ii * N_], vec1);
    dldt[ii] -= Tz_->trace_deriv(&dadt[ii * N_]) / 2;
  }
}



/// Hessian matrix of Log-likelihood.
///
/// @param[out] d2ldt length-p*p vector of the Hessian, elememts d2ldt[0:p + ii*p] is the (ii+1)-th column of Hessian matrix.
/// @param[in] z length-N vector of observation.
/// @param[in] dzdt length-N*p vector of observation derivatives, elememts dzdt[0:N + ii*N] is the derivative of z with respect to the ii-th parameter.
/// @param[in] d2zdt length-N*p*p vector of observation derivatives, elememts d2zdt[0:N + (ii*p+jj)*N] is the second derivative of z with respect to the ii-th parameter and then the jj-th parameter.
/// @param[in] acf length-N vector of auto-covariance.
/// @param[in] dadt length-N*p vector of observation derivatives, elememts dadt[0:N + ii*N] is the derivative of acf with respect to the ii-th parameter.
/// @param[in] d2adt length-N*p*p vector of observation derivatives, elememts d2adt[0:N + (ii*p+jj)*N] is the second derivative of acf with respect to the ii-th parameter and then the jj-th parameter.
inline void NormalToeplitz::hess(double* d2ldt,
				 const double* z,
				 const double* dzdt,
				 const double* d2zdt,
				 const double* acf,
				 const double* dadt,
				 const double* d2adt) {
  Tz_->setAcf(acf);
  Tz_->solve(vec1, z);
  double ans;
  std::fill(d2ldt, d2ldt + p_ * p_, 0);
  for (int ii = 0; ii < p_; ++ii) {
    for (int jj = 0; jj <= ii; ++jj) {
      Tz_->product(vec4, vec1, &dadt[jj * N_]);
      Tz_->product(vec3, vec1, &dadt[ii * N_]);
      ans = dot_prod(&d2zdt[(ii * p_ + jj) * N_], vec1);

      Tz_->solve(vec2, vec4);
      ans -= dot_prod(&dzdt[ii * N_], vec2);
      ans += dot_prod(vec3, vec2);
      Tz_->solve(vec2, vec3);
      ans -= dot_prod(&dzdt[jj * N_], vec2);
      Tz_->solve(vec2, &dzdt[jj * N_]);
      ans += dot_prod(&dzdt[ii * N_], vec2);
      ans *= 2;
			
      Tz_->product(vec2, vec1, &d2adt[(ii * p_ + jj) * N_]);
      ans -= dot_prod(vec1, vec2);
      ans += Tz_->trace_deriv(&d2adt[(ii * p_ + jj) * N_]);
      ans -= Tz_->trace_hess(&dadt[ii * N_], &dadt[jj * N_]);

      d2ldt[ii * p_ + jj] = -ans / 2;
    }
  }

  if (p_ > 1) {
    for (int ii = 0; ii < p_; ++ii) {
      for (int jj = ii + 1; jj < p_; ++jj) {
	d2ldt[ii * p_ + jj] = d2ldt[jj * p_ + ii];
      }
    }
  }

}


/// Full Gradient of the Log-Density for Auto-Differentiation Algorithms.
///
/// @param[out] dldz length-N vector of the log-density with respect to the observation z.
/// @param[out] dlda length-N vector of the log-density with respect to the auto-covariance acf.
/// @param[in] z length-N vector of observation.
/// @param[in] acf length-N vector of auto-covariance.
inline void NormalToeplitz::grad_full(double* dldz, double* dlda,
				      const double* z, const double* acf) {
  Tz_->setAcf(acf);

  // gradient with respect to z
  Tz_->solve(dldz, z);
  for (int ii = 0; ii < N_; ii++) {
    dldz[ii] = -dldz[ii];
  }

  // gradient with respect to acf
  Tz_->solve(vec1, z);	
  vec2[0] = 1;
  std::fill(vec2 + 1, vec2 + N_, 0);
  Tz_->solve(vec3, vec2);
  double tau1 = vec3[0];
  std::fill(phi, phi + N_, 0);
  phi[0] = vec1[0];
  Tz_->product(dlda, vec1, phi, vec1); // dlda = upper.toep(Vz) %*% Vz = ip
  vec4[0] = 0;
  for (int ii = 1; ii < N_; ++ii) {
    vec4[ii] = vec3[N_ - ii];
  }
  for (int ii = 0; ii < N_; ++ii) {
    vec2[ii] = (N_ - ii) * vec3[ii];
  } // vec2 = (N_:1 * tau)

  phi[0] = vec3[0];
  Tz_->product(vec1, vec2, phi, vec3); // vec1 = upper.toep(tau) %*% (N_:1 * tau) = tr

  for (int ii = 0; ii < N_; ++ii) {
    vec2[ii] = (N_ - ii) * vec4[ii];
  } // vec2 = (N_:1 * tau2)

  phi[0] = vec4[0];
  Tz_->product(vec3, vec2, phi, vec4); // vec3 = upper.toep(tau2) %*% (N_:1 * tau2)

  for (int ii = 0; ii < N_; ++ii) {
    vec1[ii] -= vec3[ii];
    vec1[ii] /= tau1;
    dlda[ii] -= vec1[ii];
  } // vec1 = (vec1 - vec3) / tau[1] = tr, dlda = ip - tr

  dlda[0] /= 2;
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
ldens += Tz->logDet() + N_ * log2Pi;
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
for (int ii = 0; ii < p_; ++ii) {
T2_->setAcf(&dadt[ii * N_]);
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
for (int ii = 0; ii < p_; ++ii) {
for (int jj = 0; jj <= ii; ++jj) {
T2_->setAcf(&dadt[jj * N_]);
T2_->multVec(vec4, vec1);
T2_->setAcf(&dadt[ii * N_]);
T2_->multVec(vec3, vec1);

ans = dot_prod(&d2zdt[(ii * p_ + jj) * N_], vec1, N_);
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

T2_->setAcf(&d2adt[(ii * p_ + jj) * N_]);
T2_->multVec(vec2, vec1);
ans -= dot_prod(vec1, vec2, N_);
ans += Tz->traceProd(&d2adt[(ii * p_ + jj) * N_]);
ans -= Tz->traceDeriv(&dadt[ii * N_], &dadt[jj * N_]);

d2ldt[ii * p_ + jj] = -ans / 2;
}
}

if (p_ > 1) {
for (int ii = 0; ii < p_; ++ii) {
for (int jj = ii + 1; jj < p_; ++jj) {
d2ldt[ii * p_ + jj] = d2ldt[jj * p_ + ii];
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

T2_->setAcf(vec1);
T2_->mult0Vec(dlda, vec1); // dlda = upper.toep(Vz) %*% Vz = ip

vec4[0] = 0;
for (int ii = 1; ii < N_; ++ii) {
vec4[ii] = vec3[N_ - ii];
} // vec4 = tau2

for (int ii = 0; ii < N_; ++ii) {
vec2[ii] = (N_ - ii) * vec3[ii];
} // vec2 = (N_:1 * tau)

T2_->setAcf(vec3);
T2_->mult0Vec(vec1, vec2); // vec1 = upper.toep(tau) %*% (N_:1 * tau) = tr

for (int ii = 0; ii < N_; ++ii) {
vec2[ii] = (N_ - ii) * vec4[ii];
} // vec2 = (N_:1 * tau2)
T2_->setAcf(vec4);
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

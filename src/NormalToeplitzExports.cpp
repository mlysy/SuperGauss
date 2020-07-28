/// @file NormalToeplitzExports.cpp
///
/// @brief `Rcpp` wapper functions to NormalToeplitz methods.

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/NormalToeplitz.h"


/// NormalToeplitz class constructor.
///
/// @brief Creates an `Rcpp::XPtr` pointer to a NormalToeplitz object and passes the pointer address to the R session.
///
/// @param N Size of NormalToeplitz random vector.
///
/// @return The `Rcpp::XPtr` pointer to the instantiated NormalToeplitz object.
//[[Rcpp::export]]
SEXP NormalToeplitz_ctor(int N) {
  NormalToeplitz *NTz = new NormalToeplitz(N);
  XPtr<NormalToeplitz> NTz_ptr(NTz, true);
  return NTz_ptr;
}

/// Log-density of NormalToeplitz distribution.
///
/// @param[in] NTz_ptr `Rcpp::XPtr` pointer to a NormalToeplitz object.
/// @param[in] z Observation matrix of size `N x n_z`, where each column is an observation.
/// @param[in] acf Autocorrelation vector of length `N`.
///
/// @return Vector value of the log-density at each column of `z`.
//[[Rcpp::export]]
NumericVector NormalToeplitz_logdens(SEXP NTz_ptr, NumericMatrix z,
				     NumericVector acf) {
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  int n_z = z.ncol();
  int N = z.nrow();
  NumericVector ldens(n_z);
  double *z_ = REAL(z);
  NTz->set_acf(REAL(acf));
  for(int ii=0; ii<n_z; ii++) {
    NTz->set_z(&z_[ii*N]);
    REAL(ldens)[ii] = NTz->logdens();
  }
  return ldens;
  // return NTz->logdens(REAL(z), REAL(acf));
}

/// Gradient of NormalToeplitz loglikelihood.
///
/// Calculates the gradient with respect to `theta` of the loglikelihood corresponding to
///
/// ```
/// z(theta) ~ NormalToeplitz(acf(theta)).
/// ```
///
/// @param[in] NTz_ptr `Rcpp::XPtr` pointer to a NormalToeplitz object.
/// @param[in] z Observation vector of length `N`.
/// @param[in] dzdt Gradient of `z` with respect to `theta`.  A matrix of size `N x n_theta`.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[in] dadt Gradient of `acf` with respect to `theta`.  A matrix of size `N x n_theta`.
/// @param[in] full_out If `true`, returns the log-density as well.
///
/// @return Gradient of the loglikelihood.  A vector of length `n_theta`, or a list with elements `ldens` for the log-density (a scalar) and `grad` for the gradient vector.
//[[Rcpp::export]]
SEXP NormalToeplitz_grad(SEXP NTz_ptr,
			 NumericVector z,
			 NumericMatrix dzdt, 
			 NumericVector acf,
			 NumericMatrix dadt,
			 bool full_out = false) {
  int n_theta = dzdt.ncol();
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  NumericVector dldt(n_theta);
  NTz->set_acf(REAL(acf));
  NTz->set_z(REAL(z));
  NTz->grad(REAL(dldt), REAL(dzdt), REAL(dadt), n_theta);
  if(full_out) {
    double ldens = NTz->logdens();
    return List::create(_["ldens"] = ldens, _["grad"] = dldt);
  } else {
    return dldt;
  }
}

/// Hessian of NormalToeplitz loglikelihood.
///
/// Calculates the Hessian matrix with respect to `theta` of the loglikelihood corresponding to
///
/// ```
/// z(theta) ~ NormalToeplitz(acf(theta)).
/// ```
///
/// @param[in] NTz_ptr `Rcpp::XPtr` pointer to a NormalToeplitz object.
/// @param[in] z Observation vector of length `N`.
/// @param[in] dzdt Gradient of `z` with respect to `theta`.  A matrix of size `N x n_theta`.
/// @param[in] d2zdt Hessian of `z` with respect to `theta`.  A matrix of size `N x (n_theta * n_theta)` corresponding to the Hessian tensor of size `N x n_theta x n_theta` flattened in column-major order.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[in] dadt Gradient of `acf` with respect to `theta`.  A matrix of size `N x n_theta`.
/// @param[in] d2adt Hessian of `acf` with respect to `theta`.  A matrix of size `N x (n_theta * n_theta)` corresponding to the Hessian tensor of size `N x n_theta x n_theta` flattened in column-major order.
/// @param[in] full_out If `true`, returns the log-density and gradient as well.///
/// @return Hessian of the loglikelihood.  A matrix of size `n_theta x n_theta`, or a list with elements `ldens`, `grad`, and `hess` consisting of the log-density (scalar), the gradient vector (of size `n_theta`) and the hessian matrix, respectively.
//[[Rcpp::export]]
SEXP NormalToeplitz_hess(SEXP NTz_ptr,
			 NumericVector z,
			 NumericMatrix dzdt,
			 NumericMatrix d2zdt,
			 NumericVector acf,
			 NumericMatrix dadt,
			 NumericMatrix d2adt,
			 bool full_out = false) {
  int n_theta = dzdt.ncol();
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  NumericMatrix d2ldt(n_theta, n_theta);
  NTz->set_acf(REAL(acf));
  NTz->set_z(REAL(z));
  NTz->hess(REAL(d2ldt), REAL(dzdt), REAL(d2zdt), 
	    REAL(dadt), REAL(d2adt), n_theta);
  if(full_out) {
    double ldens = NTz->logdens();
    NumericVector dldt(n_theta);
    NTz->grad(REAL(dldt), REAL(dzdt), REAL(dadt), n_theta);
    return List::create(_["ldens"] = ldens,
			_["grad"] = dldt,
			_["hess"] = d2ldt);
  } else {
    return d2ldt;
  }
}

/// Full gradient of NormalToeplitz log-density.
///
/// Calculates the gradient with respect to each element of `z` and `acf` of the log-density corresponding to `z ~ NormalToeplitz(acf)`.
///
/// @param[in] NTz_ptr `Rcpp::XPtr` pointer to a NormalToeplitz object.
/// @param[in] z Observation vector of length `N`.
/// @param[in] acf Autocorrelation vector of length `N`.
/// @param[in] calc_dldz Whether to calculate the gradient with respect to `z`.
/// @param[in] calc_dlda Whether or to calculate the gradient with respect to `acf`.
///
/// @return `Rcpp::List` with elements:
/// - `ldens`: The log-density evaluated at `z` and `acf`.
/// - `dldz`: The gradient with respect to `z`, if `calc_dldz = true`.
/// - `dlda`: The gradient with respect to `acf`, if `calc_dlda = true`.
//
//[[Rcpp::export]]
List NormalToeplitz_grad_full(SEXP NTz_ptr,
			      NumericVector z, NumericVector acf,
			      bool calc_dldz = true, bool calc_dlda = true) {
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  int N = NTz->size();
  double ldens;
  NumericVector dldz(calc_dldz ? N : 1);
  NumericVector dlda(calc_dlda ? N : 1);
  NTz->set_acf(REAL(acf));
  NTz->set_z(REAL(z));
  ldens = NTz->grad_full(REAL(dldz), REAL(dlda), calc_dldz, calc_dlda);
  List out;
  out["ldens"] = ldens;
  if(calc_dldz) out["dldz"] = dldz;
  if(calc_dlda) out["dlda"] = dlda;
  return out;
}

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
//[[Rcpp::export(".NormalToeplitz_constructor")]]
SEXP NormalToeplitz_constructor(int N) {
  NormalToeplitz *NTz = new NormalToeplitz(N);
  XPtr<NormalToeplitz> NTz_ptr(NTz, true);
  return NTz_ptr;
}

/// Log-density of NormalToeplitz distribution.
///
/// @param[in] NTz_ptr `Rcpp::XPtr` pointer to a NormalToeplitz object.
/// @param[in] z Observation vector of length `N`.
/// @param[in] acf Autocorrelation vector of length `N`.
///
/// @return Scalar value of the log-density.
//[[Rcpp::export(".NormalToeplitz_logdens")]]
double NormalToeplitz_logdens(SEXP NTz_ptr, NumericVector z,
			      NumericVector acf) {
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  return NTz->logdens(REAL(z), REAL(acf));
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
///
/// @return Gradient of the loglikelihood.  A vector of length `n_theta`.
//[[Rcpp::export(".NormalToeplitz_grad")]]
NumericVector NormalToeplitz_grad(SEXP NTz_ptr,
				  NumericVector z,
				  NumericMatrix dzdt, 
				  NumericVector acf,
				  NumericMatrix dadt,
				  int n_theta) {
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  NumericVector dldt(n_theta);
  NTz->grad(REAL(dldt), REAL(z), REAL(dzdt),
	    REAL(acf), REAL(dadt), n_theta);
  return dldt;
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
///
/// @return Hessian of the loglikelihood.  A matrix of size `n_theta x n_theta`.
//[[Rcpp::export(".NormalToeplitz_hess")]]
NumericMatrix NormalToeplitz_hess(SEXP NTz_ptr,
				  NumericVector z,
				  NumericMatrix dzdt,
				  NumericMatrix d2zdt,
				  NumericVector acf,
				  NumericMatrix dadt,
				  NumericMatrix d2adt,
				  int n_theta) {
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  NumericMatrix d2ldt(n_theta, n_theta);
  NTz->hess(REAL(d2ldt), REAL(z), REAL(dzdt), REAL(d2zdt), 
	    REAL(acf), REAL(dadt), REAL(d2adt), n_theta);
  return d2ldt;
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
/// @return `Rcpp::List` with one or both elements `dldz` and `dlda`, corresponding to the length `N` gradient vectors with respect to `z` and `acf`, respectively.
//
//[[Rcpp::export(".NormalToeplitz_grad_full")]]
List NormalToeplitz_grad_full(SEXP NTz_ptr,
			      NumericVector z, NumericVector acf,
			      bool calc_dldz = true, bool calc_dlda = true) {
  XPtr<NormalToeplitz> NTz(NTz_ptr);
  int N = NTz->size();
  NumericVector dldz(calc_dldz ? N : 1);
  NumericVector dlda(calc_dlda ? N : 1);
  NTz->grad_full(REAL(dldz), REAL(dlda), REAL(z), REAL(acf),
		 calc_dldz, calc_dlda);
  List out;
  if(calc_dldz) out["dldz"] = dldz;
  if(calc_dlda) out["dlda"] = dlda;
  // List::create(Named("dldz") = dldz, _["dlda"] = dlda);
  return out;
}

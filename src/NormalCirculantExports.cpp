/// @file NormalCirculantExports.cpp
///
/// @brief `Rcpp` wapper functions to NormalCirculant methods.

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/NormalCirculant.h"


/// NormalCirculant class constructor.
///
/// @brief Creates an `Rcpp::XPtr` pointer to a NormalCirculant object and passes the pointer address to the R session.
///
/// @param N Size of NormalCirculant random vector.
///
/// @return The `Rcpp::XPtr` pointer to the instantiated NormalCirculant object.
// [[Rcpp::export]]
SEXP NormalCirculant_ctor(int N) {
  NormalCirculant *NCt = new NormalCirculant(N);
  XPtr<NormalCirculant> pNCt(NCt, true);
  return pNCt;
}

/// Log-density of NormalCirculant distribution.
///
/// @param[in] pNCt `Rcpp::XPtr` pointer to a NormalCirculant object.
/// @param[in] z Observation matrix of size `N x n_z`, where each column is an observation.
/// @param[in] uacf Half-autocorrelation vector of length `Nu = floor(N/2)+1`.
///
/// @return Vector value of the log-density at each column of `z`.
// [[Rcpp::export]]
NumericVector NormalCirculant_logdens(SEXP pNCt, NumericMatrix z,
				      NumericVector uacf) {
  XPtr<NormalCirculant> NCt(pNCt);
  int n_z = z.ncol();
  int N = z.nrow();
  NumericVector ldens(n_z);
  double *z_ = REAL(z);
  NCt->set_acf(REAL(uacf));
  for(int ii=0; ii<n_z; ii++) {
    NCt->set_z(&z_[ii*N]);
    REAL(ldens)[ii] = NCt->logdens();
  }
  return ldens;
  // return NCt->logdens(REAL(z), REAL(uacf));
}

/// Full gradient of NormalCirculant log-density.
///
/// Calculates the gradient with respect to each element of `z` and `uacf` of the log-density corresponding to `z ~ NormalCirculant(uacf)`.
///
/// @param[in] pNCt `Rcpp::XPtr` pointer to a NormalCirculant object.
/// @param[in] z Observation vector of length `N`.
/// @param[in] uacf Half-autocorrelation vector of length `Nu = floor(N/2)+1`.
/// @param[in] calc_dldz Whether to calculate the gradient with respect to `z`.
/// @param[in] calc_dldu Whether or to calculate the gradient with respect to `uacf`.
///
/// @return `Rcpp::List` with elements:
/// - `ldens`: The log-density evaluated at `z` and `uacf`.
/// - `dldz`: The gradient with respect to `z`, if `calc_dldz = true`.
/// - `dldu`: The gradient with respect to `uacf`, if `calc_dldu = true`.
//
// [[Rcpp::export]]
List NormalCirculant_grad_full(SEXP pNCt,
			       NumericVector z, NumericVector uacf,
			       bool calc_dldz = true,
			       bool calc_dldu = true) {
  XPtr<NormalCirculant> NCt(pNCt);
  int N = NCt->size();
  int Nu = N/2 + 1;
  double ldens;
  NumericVector dldz(calc_dldz ? N : 1);
  NumericVector dldu(calc_dldu ? Nu : 1);
  ldens = NCt->grad_full(REAL(dldz), REAL(dldu), REAL(z), REAL(uacf),
				calc_dldz, calc_dldu);
  List out;
  out["ldens"] = ldens;
  if(calc_dldz) out["dldz"] = dldz;
  if(calc_dldu) out["dldu"] = dldu;
  return out;
}

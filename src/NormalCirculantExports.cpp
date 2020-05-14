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
/// @param[in] z Observation vector of length `N`.
/// @param[in] uacf Half-autocorrelation vector of length `Nu = floor(N/2)+1`.
///
/// @return Scalar value of the log-density.
// [[Rcpp::export]]
double NormalCirculant_logdens(SEXP pNCt, NumericVector z,
			      NumericVector uacf) {
  XPtr<NormalCirculant> NCt(pNCt);
  return NCt->logdens(REAL(z), REAL(uacf));
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
/// @return `Rcpp::List` with one or both elements `dldz` and `dldu`, corresponding to the length `N` gradient vector with respect to `z` and the length `Nu = floor(N/2)+1` gradient vector with respect to `uacf`.
//
// [[Rcpp::export]]
List NormalCirculant_grad_full(SEXP pNCt,
			       NumericVector z, NumericVector uacf,
			       bool calc_dldz = true,
			       bool calc_dldu = true) {
  XPtr<NormalCirculant> NCt(pNCt);
  int N = NCt->size();
  int Nu = N/2 + 1;
  NumericVector dldz(calc_dldz ? N : 1);
  NumericVector dldu(calc_dldu ? Nu : 1);
  NCt->grad_full(REAL(dldz), REAL(dldu), REAL(z), REAL(uacf),
		 calc_dldz, calc_dldu);
  List out;
  if(calc_dldz) out["dldz"] = dldz;
  if(calc_dldu) out["dldu"] = dldu;
  return out;
}

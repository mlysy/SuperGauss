/// @file CirculantExports.cpp
///
/// @brief Rcpp wrappers for Circulant matrix class.

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/Circulant.h"

/// Construct a Circulant matrix object.
///
/// Instantiates a `Circulant` object on the C++ side, wraps it in an `Rcpp::XPtr`, and returns the corresponding `externalptr` on the R side.
///
/// @param[in] N Size of Circulant matrix.
/// @return An `externalptr` pointing to the Circulant object.
///
// [[Rcpp::export]]
SEXP Circulant_ctor(int N) {
  Circulant *Circ = new Circulant(N);
  XPtr<Circulant> pCirc(Circ, true);
  return pCirc;
}

/// Set the autocorrelation of the Circulant matrix.
///
/// @param[in] pCirc `externalptr` pointer to Circulant matrix. 
/// @param[in] uacf The unique elements of the autocorrelation vector.
///
// [[Rcpp::export]]
void Circulant_set_acf(SEXP pCirc, NumericVector uacf) {
  XPtr<Circulant> Circ(pCirc);
  Circ->set_acf(REAL(uacf));
  return;
}

/// Get the autocorrelation of the Circulant matrix.
///
/// @param[in] pCirc `externalptr` pointer to Circulant matrix.
/// @return The complete autocorrelation vector of length `N`.
///
// [[Rcpp::export]]
NumericVector Circulant_get_acf(SEXP pCirc) {
  XPtr<Circulant> Circ(pCirc);
  NumericVector acf(Circ->size());
  Circ->get_acf(REAL(acf));
  return acf;
}

/// Set the PSD of the Circulant matrix.
///
/// @param[in] pCirc `externalptr` pointer to Circulant matrix. 
/// @param[in] upsd The unique elements of the PSD vector.
///
// [[Rcpp::export]]
void Circulant_set_psd(SEXP pCirc, NumericVector upsd) {
  XPtr<Circulant> Circ(pCirc);
  Circ->set_psd(REAL(upsd));
  return;
}

/// Get the PSD of the Circulant matrix.
///
/// @param[in] pCirc `externalptr` pointer to Circulant matrix.
/// @return The complete PSD vector of length `N`.
///
// [[Rcpp::export]]
NumericVector Circulant_get_psd(SEXP pCirc) {
  XPtr<Circulant> Circ(pCirc);
  NumericVector psd(Circ->size());
  Circ->get_psd(REAL(psd));
  return psd;
}


/// Checks whether the acf of the Circulant matrix has been set.
///
/// @param[in] pCirc `externalptr` pointer to Circulant matrix.
/// @return Logical; whether its acf has been set.
// [[Rcpp::export]]
bool Circulant_has_acf(SEXP pCirc) {
  XPtr<Circulant> Circ(pCirc);
  return Circ->has_acf();
}

/// Circulant matrix-vector product.
///
/// @param[in] pCirc `externalptr` pointer to a Circulant matrix of size `N`.
/// @param[in] X Matrix of size `N x p`.
/// @return Output matrix of size `N x p` for the matrix-vector multiplication `Y = Circulant(acf) * X`.
// [[Rcpp::export]]
NumericVector Circulant_prod(SEXP pCirc, NumericMatrix X) {
   XPtr<Circulant> Circ(pCirc);
   int N = X.nrow();
   int p = X.ncol();
   NumericMatrix Y(N, p);
   for(int ii=0; ii<p; ii++) {
     Circ->prod(&REAL(Y)[N*ii], &REAL(X)[N*ii]);
   }
   return Y;
}


/// Solve Circulant system of equations.
///
/// @param[in] pCirc `externalptr` pointer to a Circulant matrix of size `N`.
/// @param[in] X Matrix of size `N x p`.
/// @return Output matrix of size `N x p` for the calculation of `Y = Circulant(acf)^{-1} * X`.
// [[Rcpp::export]]
NumericVector Circulant_solve(SEXP pCirc, NumericMatrix X) {
   XPtr<Circulant> Circ(pCirc);
   int N = X.nrow();
   int p = X.ncol();
   NumericMatrix Y(N, p);
   for(int ii=0; ii<p; ii++) {
     Circ->solve(&REAL(Y)[N*ii], &REAL(X)[N*ii]);
   }
   return Y;
}

/// Log-determinant of the Circulant matrix.
///
/// @param[in] pCirc `externalptr` pointer to a Circulant matrix of size `N`.
/// @return The log-determinant `log(det(Circulant(acf)))`.
// [[Rcpp::export]]
double Circulant_log_det(SEXP pCirc) {
   XPtr<Circulant> Circ(pCirc);
   return Circ->log_det();
}


/// @file ToeplitzExports.cpp
///
/// @brief Rcpp wrappers for Toeplitz matrix class.

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/Toeplitz.h"

/// Construct a Toeplitz matrix object.
///
/// Instantiates a `Toeplitz` object on the C++ side, wraps it in an `Rcpp::XPtr`, and returns the corresponding `externalptr` on the R side.
///
/// @param[in] N Size of Toeplitz matrix.
/// @return An `externalptr` pointing to the Toeplitz object.
///
// [[Rcpp::export]]
SEXP Toeplitz_ctor(int N) {
  Toeplitz *Toep = new Toeplitz(N);
  XPtr<Toeplitz> pToep(Toep, true);
  return pToep;
}

/// Set the autocorrelation of the Toeplitz matrix.
///
/// @param[in] pToep `externalptr` pointer to Toeplitz matrix. 
/// @param[in] acf Autocorrelation vector of length `N`.
///
// [[Rcpp::export]]
void Toeplitz_set_acf(SEXP pToep, NumericVector acf) {
  XPtr<Toeplitz> Toep(pToep);
  Toep->set_acf(REAL(acf));
  return;
}

/// Get the autocorrelation of the Toeoplitz matrix.
///
/// @param[in] pToep `externalptr` pointer to Toeplitz matrix.
/// @return The autocorrelation vector of length `N`.
///
// [[Rcpp::export]]
NumericVector Toeplitz_get_acf(SEXP pToep) {
  XPtr<Toeplitz> Toep(pToep);
  NumericVector acf(Toep->size());
  Toep->get_acf(REAL(acf));
  return acf;
}

/// Toeplitz matrix product.
///
/// @param[in] pToep `externalptr` pointer to a Toeplitz matrix of size `N`.
/// @param[in] X Matrix of size `N x p`.
/// @return Output matrix of size `N x p` for the matrix multiplication `Y = Toeplitz(acf) * X`.
// [[Rcpp::export]]
NumericMatrix Toeplitz_prod(SEXP pToep, NumericMatrix X) {
  XPtr<Toeplitz> Toep(pToep);
  int p = X.ncol();
  int N = X.nrow();
  NumericMatrix Y(N,p);
  for(int ii=0; ii<p; ii++) {
    Toep->prod(&REAL(Y)[N*ii], &REAL(X)[N*ii]);
  }
  return Y;
}

/// Solve Toeplitz system of equations.
///
/// @param[in] pToep `externalptr` pointer to a Toeplitz matrix of size `N`.
/// @param[in] X Matrix of size `N x p`.
/// @return Output matrix of size `N x p` for the calculation of `Y = Toeplitz(acf)^{-1} * X`.
// [[Rcpp::export]]
NumericMatrix Toeplitz_solve(SEXP pToep, NumericMatrix X) {
  XPtr<Toeplitz> Toep(pToep);
  int p = X.ncol();
  int N = X.nrow();
  NumericMatrix Y(N,p);
  for(int ii=0; ii<p; ii++) {
    Toep->solve(&REAL(Y)[N*ii], &REAL(X)[N*ii]);
  }
  return Y;
}

/// Log-determinant of the Toeplitz matrix.
///
/// @param[in] pToep `externalptr` pointer to a Toeplitz matrix of size `N`.
/// @return The log-determinant `log(det(Toeplitz(acf)))`.
// [[Rcpp::export]]
double Toeplitz_log_det(SEXP pToep) {
  XPtr<Toeplitz> Toep(pToep);
  return Toep->log_det();
}

/// Gradient-specialized trace-product.
///
/// Computes `trace( Tz^{-1} * Tz0 )`, where `Tz = Toeplitz(acf)` and `Tz0 = Toeplitz(acf0)`.  This trace-product appears in the computation of
/// ```
///  d/dx log(det(Tz)) = trace( Tz^{-1} * Toeplitz(d/dx Tz),
/// ```
/// i.e., where `Tz = Toeplitz(acf(x))` is a function of `x`.
///
/// @param[in] pToep `externalptr` pointer to a Toeplitz matrix of size `N`.
/// @param[in] acf0 Vector of length `N` giving the first row/column of the Toeplitz matrix `Tz0 = Toeplitz(acf0)`.
/// @return The Toeplitz trace-product `trace( Tz^{-1} * Tz0 )`.
///
// [[Rcpp::export]]
double Toeplitz_trace_grad(SEXP pToep, NumericVector acf0) {
  XPtr<Toeplitz> Toep(pToep);
  return Toep->trace_grad(REAL(acf0));
}

/// Hessian-specialized trace-product.
///
/// Computes `trace((Tz^{-1} * Tz1) * (Tz^{-1} * Tz2))`, where `Tz = Toeplitz(acf)`, `Tz1 = Toeplitz(acf1)`, and `Tz2 = Toeplitz(acf2)`.  This trace-product appears in the computation of
/// ```
/// d^2/dxdy log(det(Tz)) = trace(Tz^{-1} d^2/dxdy Tz) - trace((Tz^{-1} d/dx Tz) * (Tz^{-1} d/dy Tz)),
/// ```
/// i.e., where `Tz = Toeplitz(acf(x,y))` is a function of `x` and `y`.
///
/// @param[in] pToep `externalptr` pointer to a Toeplitz matrix of size `N`.
/// @param[in] acf1 Vector of length `N` giving the first row/column of the  Toeplitz matrix `Tz1 = Toeplitz(acf1)`.
/// @param[in] acf2 Vector of length `N` giving the first row/column of the  Toeplitz matrix `Tz2 = Toeplitz(acf2)`.
/// @return The Toeplitz trace-product `trace((Tz^{-1} * Tz1) * (Tz^{-1} * Tz2))`.
///
// [[Rcpp::export]]
double Toeplitz_trace_hess(SEXP pToep,
			   NumericVector acf1, NumericVector acf2) {
  XPtr<Toeplitz> Toep(pToep);
  return Toep->trace_hess(REAL(acf1), REAL(acf2));
}

/// Checks whether the acf of the Toeplitz matrix has been set.
///
/// @param[in] pToep `externalptr` pointer to Toeplitz matrix.
/// @return Logical; whether its acf has been set.
// [[Rcpp::export]]
bool Toeplitz_has_acf(SEXP pToep) {
  XPtr<Toeplitz> Toep(pToep);
  return Toep->has_acf();
}


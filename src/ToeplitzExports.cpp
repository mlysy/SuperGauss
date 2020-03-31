#include <Rcpp.h>
using namespace Rcpp;
#include "Toeplitz.h"

// R wapper functions to Toeplitz methods using XPtr

//[[Rcpp::export(".Toeplitz_constructor")]]
SEXP Toeplitz_constructor(int n) {
  Toeplitz *Toep = new Toeplitz(n);
  XPtr<Toeplitz> Toep_ptr(Toep, true);
  return Toep_ptr;
}

//[[Rcpp::export(".Toeplitz_set_acf")]]
void Toeplitz_set_acf(SEXP Toep_ptr, NumericVector acf) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  Toep->set_acf(REAL(acf));
  return;
}

//[[Rcpp::export(".Toeplitz_get_acf")]]
NumericVector Toeplitz_get_acf(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  NumericVector acf(Toep->size());
  Toep->get_acf(REAL(acf));
  return acf;
}

// //[[Rcpp::export(".Toeplitz_getPhi")]]
// NumericVector Toeplitz_getPhi(SEXP Toep_ptr) {
//   XPtr<Toeplitz> Toep(Toep_ptr);
//   NumericVector phi(Toep->size());
//   Toep->getPhi(REAL(phi));
//   return phi;
// }

//[[Rcpp::export(".Toeplitz_Multiply")]]
NumericMatrix Toeplitz_Multiply(SEXP Toep_ptr, NumericMatrix X) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  int p = X.ncol();
  int n = X.nrow();
  NumericMatrix Y(n,p);
  for(int ii=0; ii<p; ii++) {
    Toep->mult(&REAL(Y)[n*ii], &REAL(X)[n*ii]);
  }
  return Y;
}

//[[Rcpp::export(".Toeplitz_Solve")]]
NumericMatrix Toeplitz_Solve(SEXP Toep_ptr, NumericMatrix X) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  int p = X.ncol();
  int n = X.nrow();
  NumericMatrix Y(n,p);
  for(int ii=0; ii<p; ii++) {
    Toep->solve(&REAL(Y)[n*ii], &REAL(X)[n*ii]);
  }
  return Y;
}

//[[Rcpp::export(".Toeplitz_log_det")]]
double Toeplitz_log_det(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  return Toep->log_det();
}

// //[[Rcpp::export(".Toeplitz_traceT2")]]
// double Toeplitz_traceT2(SEXP Toep_ptr, NumericVector acf2) {
//   XPtr<Toeplitz> Toep(Toep_ptr);
//   return Toep->traceProd(REAL(acf2));
// }

// //[[Rcpp::export(".Toeplitz_traceT4")]]
// double Toeplitz_traceT4(SEXP Toep_ptr, NumericVector acf2, NumericVector acf3) {
//   XPtr<Toeplitz> Toep(Toep_ptr);
//   return Toep->traceDeriv(REAL(acf2), REAL(acf3));
// }

//[[Rcpp::export(".Toeplitz_has_acf")]]
bool Toeplitz_has_acf(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  return Toep->has_acf();
}


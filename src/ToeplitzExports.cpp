#include <Rcpp.h>
using namespace Rcpp;
#include "Toeplitz.h"

// R wapper functions to Toeplitz methods using XPtr

//[[Rcpp::export(".Toeplitz_constructor")]]
SEXP Toeplitz_constructor(int n, int b) {
  Toeplitz *Toep = new Toeplitz(n, b);
  XPtr<Toeplitz> Toep_ptr(Toep, true);
  return Toep_ptr;
}

//[[Rcpp::export(".Toeplitz_setAcf")]]
void Toeplitz_setAcf(SEXP Toep_ptr, NumericVector acf) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  Toep->setAcf(REAL(acf));
  return;
}

//[[Rcpp::export(".Toeplitz_getAcf")]]
NumericVector Toeplitz_getAcf(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  NumericVector acf(Toep->size());
  Toep->getAcf(REAL(acf));
  return acf;
}

//[[Rcpp::export(".Toeplitz_getPhi")]]
NumericVector Toeplitz_getPhi(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  NumericVector phi(Toep->size());
  Toep->getPhi(REAL(phi));
  return phi;
}

//[[Rcpp::export(".Toeplitz_Multiply")]]
NumericMatrix Toeplitz_Multiply(SEXP Toep_ptr, NumericMatrix X) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  int p = X.ncol();
  int n = X.nrow();
  NumericMatrix Y(n,p);
  for(int ii=0; ii<p; ii++) {
    Toep->multVec(&REAL(Y)[n*ii], &REAL(X)[n*ii]);
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
    Toep->solveVec(&REAL(Y)[n*ii], &REAL(X)[n*ii]);
  }
  return Y;
}

//[[Rcpp::export(".Toeplitz_Determinant")]]
double Toeplitz_Determinant(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  return Toep->logDet();
}

//[[Rcpp::export(".Toeplitz_traceT2")]]
double Toeplitz_traceT2(SEXP Toep_ptr, NumericVector acf2) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  return Toep->traceProd(REAL(acf2));
}

//[[Rcpp::export(".Toeplitz_traceT4")]]
double Toeplitz_traceT4(SEXP Toep_ptr, NumericVector acf2, NumericVector acf3) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  return Toep->traceDeriv(REAL(acf2), REAL(acf3));
}

//[[Rcpp::export(".Toeplitz_hasAcf")]]
bool Toeplitz_hasAcf(SEXP Toep_ptr) {
  XPtr<Toeplitz> Toep(Toep_ptr);
  return Toep->hasAcf();
}


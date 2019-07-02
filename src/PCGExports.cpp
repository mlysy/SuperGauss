#include <Rcpp.h>
using namespace Rcpp;
#include "PCG.h"

// R wapper functions to Toeplitz methods using XPtr

//[[Rcpp::export(".PCG_constructor")]]
SEXP PCG_constructor(int n) {
  PCG *P1 = new PCG(n);
  XPtr<PCG> PCG_ptr(P1, true);
  return PCG_ptr;
}

//[[Rcpp::export(".PCG_solve")]]
NumericVector PCG_Solve(SEXP PCG_ptr, NumericVector acf, NumericVector y, double tol) {
  XPtr<PCG> P1(PCG_ptr);
  int n = y.size();
  NumericVector Y(n);
  P1->solve(REAL(Y), REAL(acf), REAL(y), tol);
  return Y;
}

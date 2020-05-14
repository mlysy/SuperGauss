#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/PCG.h"

// R wapper functions to Toeplitz methods using XPtr

//[[Rcpp::export(".PCG_constructor")]]
SEXP PCG_constructor(int n) {
  PCG *P1 = new PCG(n);
  XPtr<PCG> PCG_ptr(P1, true);
  return PCG_ptr;
}

//[[Rcpp::export(".PCG_solve")]]
NumericMatrix PCG_solve(SEXP PCG_ptr, NumericVector acf, NumericMatrix y, double tol) {
  XPtr<PCG> P1(PCG_ptr);
  int N = y.nrow();
  int p = y.ncol();
  NumericMatrix Y(N, p);
  for(int ii=0; ii<p; ii++) {
    P1->solve(REAL(Y)+ii*N, REAL(acf), REAL(y)+ii*N, tol);
  }
  return Y;
}

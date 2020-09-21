/// @file PCGExports.cpp
///
/// @brief `Rcpp` wrapper functions to PCG methods.

#include <Rcpp.h>
using namespace Rcpp;
#include "SuperGauss/PCG.h"

/// PCG class constructor.
///
/// @brief Creates an `Rcpp::XPtr` pointer to a PCG object and passes the pointer address to the R session.
///
/// @param N Size of Toeplitz matrix.
///
/// @return The `Rcpp::XPtr` pointer to the instantiated PCG object.
// [[Rcpp::export]]
SEXP PCG_ctor(int N) {
  PCG *PTz = new PCG(N);
  XPtr<PCG> pPTz(PTz, true);
  return pPTz;
}

/// Solve Toeplitz system of equations using the preconditioned conjugate gradient (PCG) method.
///
/// Calculates `Y = Tz^{-1} X`, where `Tz = Toeplitz(acf)`.
///
/// @param[in] pPTz `Rcpp::XPtr` pointer to a PCG object.
/// @param[in] acf Vector of length `N` defining the Toeplitz matrix.
/// @param[in] x Matrix of size `N x p` defining the RHS of the system.
/// @param[in] tol Positive scalar specifying the tolerance of the conjugate gradient algorithm.
///
/// @return Matrix of size `N x p` containing the solution to the system.
///
// [[Rcpp::export]]
NumericMatrix PCG_solve(SEXP pPTz, NumericVector acf, NumericMatrix y, double tol) {
  XPtr<PCG> PTz(pPTz);
  int N = y.nrow();
  int p = y.ncol();
  NumericMatrix Y(N, p);
  for(int ii=0; ii<p; ii++) {
    PTz->solve(REAL(Y)+ii*N, REAL(acf), REAL(y)+ii*N, tol);
  }
  return Y;
}

/// @file mem_test.cpp
///
/// @brief R6 memory allocation tests.

// [[Rcpp::depends(SuperGauss)]]
#include <VectorFFT.h>
#include <GSchur.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP VectorFFT_ctor(int N) {
  VectorFFT *vf = new VectorFFT(N);
  XPtr<VectorFFT> pvf(vf, true);
  return pvf;
}

// [[Rcpp::export]]
SEXP GSchur2K_ctor(int N) {
  GSchur2K *g2k = new GSchur2K(N);
  XPtr<GSchur2K> pg2k(g2k, true);
  return pg2k;  
}

// [[Rcpp::export]]
SEXP GSchurN_ctor(int N) {
  GSchurN *gn = new GSchurN(N);
  XPtr<GSchurN> pgn(gn, true);
  return pgn;  
}

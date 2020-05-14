// Durbin Levinson algorithms
//
// (1) inner products of the form:
//     M = X' toeplitz(acf)^{-1} Y
// (2) obtain either X or Z from the other in:
//     X = chol(toeplitz(acf))' Z

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "SuperGauss/DurbinLevinsonXZ.h"
#include "SuperGauss/DurbinLevinsonIP.h"

////////////////////////////////////////////////////

// DurbinLevinsonXZ

// different functions for X->Z and Z->X

// toeplitzXZ
//[[Rcpp::export]]
Eigen::MatrixXd DurbinLevinson_XZ(Eigen::MatrixXd X, Eigen::VectorXd acf) {
  int n, k;
  n = acf.size();
  k = X.cols();
  // output
  MatrixXd Z(n,k);
  // tmp
  VectorXd phi(n);
  VectorXd phi2(n);
  VectorXd res(k);
  MatrixXd Xt(k,n);
  MatrixXd Zt(k,n);

  Xt = X.transpose();
  DurbinLevinsonXZ(Xt, Zt, acf, phi, phi2, res, false);
  Z = Zt.transpose();
  return(Z);
}

// toeplitzZX
//[[Rcpp::export]]
Eigen::MatrixXd DurbinLevinson_ZX(Eigen::MatrixXd Z, Eigen::VectorXd acf) {
  int n, k;
  n = acf.size();
  k = Z.cols();
  // output
  MatrixXd X(n,k);
  // tmp
  VectorXd phi(n);
  VectorXd phi2(n);
  VectorXd res(k);
  MatrixXd Xt(k,n);
  MatrixXd Zt(k,n);

  Zt = Z.transpose();
  DurbinLevinsonXZ(Xt, Zt, acf, phi, phi2, res, true);
  X = Xt.transpose();
  return(X);
}

////////////////////////////////////////////////////

// DurbinLevinsonIP

// Two versions: Eigen and Base.  Former is almost always faster.

// DurbinLevinsonEigen
//[[Rcpp::export]]
Rcpp::List DurbinLevinson_Eigen(Eigen::MatrixXd X, Eigen::MatrixXd Y,
				     Eigen::VectorXd acf,
				     int calcMode = 1) {
  int n, d, k;
  n = acf.size();
  d = X.cols();
  k = calcMode == 1 ? d : Y.cols();
  // output
  MatrixXd M(calcMode != 2 ? d : 1, k);
  double ldV = 0.0;
  // tmp
  VectorXd phi(n);
  VectorXd phi2(n);
  VectorXd rx(d);
  VectorXd ry(k);
  MatrixXd Xt(d,n);
  MatrixXd Yt(calcMode == 1 ? 1 : k, calcMode == 1 ? 1 : n);

  Xt = X.transpose();
  if(calcMode != 1) {
    Yt = Y.transpose();
  }
  DurbinLevisonEigen(M, ldV, Xt, Yt, acf, phi, phi2, rx, ry, calcMode);
  return List::create(_["IP"] = wrap(M), _["ldV"] = wrap(ldV));
}

// DurbinLevinsonBase
//[[Rcpp::export]]
Rcpp::List DurbinLevinson_Base(NumericMatrix X, NumericMatrix Y,
			       NumericVector acf, int calcMode = 1) {
  int n, d, k;
  n = acf.length();
  d = X.ncol();
  k = calcMode == 1 ? d : Y.ncol();
  // output
  NumericMatrix M(calcMode != 2 ? d : 1, k);
  double ldV = 0.0;
  // tmp
  double *phi = new double[n];
  double *phi2 = new double[n];
  double *rx = new double[d];
  double *ry = new double[k];

  DurbinLevinsonBase(REAL(M), ldV, REAL(X), REAL(Y), REAL(acf),
		     phi, phi2, rx, ry, n, d, k, calcMode);
  delete [] phi;
  delete [] phi2;
  delete [] rx;
  delete [] ry;
  return List::create(_["IP"] = wrap(M), _["ldV"] = wrap(ldV));
}

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
#include "DurbinLevinsonXZ.h"

////////////////////////////////////////////////////

// Eigen version
// NOTE: X and Y are passed as d x n and k x n matrices
template <typename T1>
void DurbinLevisonEigen(const Eigen::MatrixBase<T1>& M, double &ldV,
        const Ref <const MatrixXd>& X, const Ref <const MatrixXd>& Y,
        const Ref <const VectorXd>& acf,
        Ref <VectorXd> phi, Ref <VectorXd> phi2,
        Ref <VectorXd> rx, Ref <VectorXd> ry,
        int calcMode) {
#define _M const_cast<Eigen::MatrixBase<T1>& >(M)
  int n = acf.size();
  double nu, rp;
  // initialize
  nu = acf(0);
  ldV = 0.0;
  _M.setZero();
  for(int ii = 0; ii < n; ii++) {
    // variance
    if(ii > 0) nu *= (1.0-phi(ii-1)*phi(ii-1));
    ldV += log(nu);
    // residuals
    rx = X.leftCols(ii) * phi.head(ii).reverse();
    rx = X.col(ii) - rx;
    if(calcMode != 1) {
      ry = Y.leftCols(ii) * phi.head(ii).reverse();
      ry = Y.col(ii) - ry;
    }
    // inner products
    if(calcMode == 0) {
      // matrix inner product (default)
      rx /= nu;
      _M.noalias() += rx * ry.transpose();
    }
    else if(calcMode == 1) {
      // square norm (i.e. lower triangular only)
      rx /= sqrt(nu);
      _M.template selfadjointView<Eigen::Lower>().rankUpdate(rx);
    } else {
      // elementwise inner product (i.e. diag(<X, Y>) )
      _M.array() += (rx.array() * ry.array())/nu;
    }
    // new coefficients
    if(ii < n-1) {
      // coefficients
      phi2.head(ii) = phi.head(ii).reverse();
      rp = phi2.head(ii).dot(acf.segment(1,ii));
      //rp = phi.head(ii).dot(acf.segment(1,ii).reverse());
      //phi2.head(ii) = phi.head(ii).reverse();
      phi(ii) = (acf(ii+1) - rp)/nu;
      phi.head(ii) -= phi(ii) * phi2.head(ii);
    }
  }
  if(calcMode == 1) {
    // upper triangular part of M
    _M.template triangularView<Eigen::Upper>() = M.transpose();
  }
  //std::cout << M << std::endl;
#undef _M
  return;
}

// Base version
// NOTE: X and Y are passed as n x d and n x k matrices,
// but in columnwise pure double format
void DurbinLevinsonBase(double *M, double &ldV, double *X,
			double *Y, double *acf,
			double *phi, double *phi2, double *rx,
			double *ry, int n, int d, int k, int calcMode) {
  int xI, yJ, ii, jj;
  double nu, rp;
  nu = acf[0];
  ldV = 0.0;
  for(ii = 0; ii < n; ii++) {
    // variance
    if(ii > 0) nu *= (1-phi[ii-1]*phi[ii-1]);
    ldV += log(nu);
    // residuals
    for(xI = 0; xI < d; xI++) rx[xI] = 0.0;
    for(yJ = 0; yJ < k; yJ++) ry[yJ] = 0.0;
    for(jj = 0; jj < ii; jj++) {
      for(xI = 0; xI < d; xI++) rx[xI] += phi[jj] * X[xI*n + ii-1-jj];
      for(yJ = 0; yJ < k; yJ++) ry[yJ] += phi[jj] * Y[yJ*n + ii-1-jj];
    }
    for(xI = 0; xI < d; xI++) rx[xI] = X[xI*n + ii] - rx[xI];
    for(yJ = 0; yJ < k; yJ++) ry[yJ] = Y[yJ*n + ii] - ry[yJ];
    // inner products
    if(calcMode == 0) {
      // matrix inner product (default)
      for(yJ = 0; yJ < k; yJ++) {
        for(xI = 0; xI < d; xI++) M[yJ*d + xI] += rx[xI]*ry[yJ]/nu;
      }
    }
    else if(calcMode == 1) {
      // square norm (i.e. lower triangular only)
      for(yJ = 0; yJ < d; yJ++) {
        M[yJ*d + yJ] += rx[yJ]*rx[yJ]/nu;
        for(xI = 0; xI < yJ; xI++) {
          M[yJ*d + xI] += rx[xI]*ry[yJ]/nu;
        }
      }
    } else {
      // elementwise inner product (i.e. diag(<X, Y>) )
      for(xI = 0; xI < d; xI++) M[xI] += rx[xI]*ry[xI]/nu;
    }
    // new coefficients
    if(ii < n-1) {
      // coefficients
      rp = 0.0;
      for(jj = 0; jj < ii; jj++) {
        rp += phi[jj] * acf[ii-jj];
        phi2[jj] = phi[ii-1-jj];
      }
      phi[ii] = (acf[ii+1] - rp)/nu;
      for(jj = 0; jj < ii; jj++) phi[jj] -= phi[ii]*phi2[jj];
    }
  }
  if(calcMode == 1) {
    for(yJ = 0; yJ < d; yJ++) {
      for(xI = 0; xI < yJ; xI++) M[xI*d + yJ] = M[yJ*d + xI];
    }
  }
}

// DurbinLevinsonXZ

// different functions for X->Z and Z->X

//[[Rcpp::export("toeplitzXZ")]]
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

//[[Rcpp::export("toeplitzZX")]]
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

//[[Rcpp::export("DurbinLevinsonEigen")]]
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

//[[Rcpp::export("DurbinLevinsonBase")]]
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

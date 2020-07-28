// inner product calculation by Durbin-Levinson algorithm
//
// calcMode = 0 for M = X' iT Y, where iT = toeplitz(acf)^{-1}
// calcMode = 1 for M = X' iT X,
// calcMode = 2 for M = diag(X' iT X)

#ifndef DurbinLevinsonIP_h
#define DurbinLevinsonIP_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// Eigen version
// NOTE: X and Y are passed as d x n and k x n matrices
template <typename T1>
inline void DurbinLevisonEigen(const Eigen::MatrixBase<T1>& M, double &ldV,
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
inline void DurbinLevinsonBase(double *M, double &ldV, double *X,
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

#endif

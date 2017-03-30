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
void DurbinLevinsonEigen(const Eigen::MatrixBase<T1>& M, double &ldV,
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
void DurbinLevinsonBase(double *M, double &ldV, double *X,
			double *Y, double *acf,
			double *phi, double *phi2, double *rx,
			double *ry, int n, int d, int k, int calcMode);


#endif

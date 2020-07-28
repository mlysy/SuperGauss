// conversion between X and Z by Durbin-Levinson algorithm
//
// Here, X = chol(toeplitz(acf))' Z,
// and we can solve for either X or Z having supplied the other using
// the boolean argument ZtoX

#ifndef DurbinLevinsonXZ_h
#define DurbinLevinsonXZ_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


// NOTE: inputs are transposed versions of algebra X and Z.  That is,
// X and Z are m x n matrices, where n = length(acf).
template <typename T1, typename T2>
inline void DurbinLevinsonXZ(const Eigen::MatrixBase<T1>& X,
        const Eigen::MatrixBase<T2>& Z,
        const Ref <const VectorXd>& acf,
        Ref <VectorXd> phi, Ref <VectorXd> phi2, Ref <VectorXd> res,
        bool ZtoX) {
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _Z const_cast<Eigen::MatrixBase<T2>& >(Z)
  int n; //,m;
  n = acf.size();
  //m = X.rows();
  double nu, sqrtNu, rp;
  nu = acf(0);
  int ii; //, jj, kk;
  for(ii = 0; ii < n; ii++) {
    // variance
    if(ii > 0) nu *= (1-phi(ii-1)*phi(ii-1));
    sqrtNu = sqrt(nu);
    // mean and actual conversion
    if(ii == 0) {
      res.setZero();
    }
    else {
      res = _X.leftCols(ii) * phi.head(ii).reverse();
    }
    if(ZtoX) {
      _X.col(ii) = res + sqrtNu * _Z.col(ii);
    }
    else {
      _Z.col(ii) = (_X.col(ii) - res)/sqrtNu;
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
#undef _X
#undef _Z
  return;
}

#endif

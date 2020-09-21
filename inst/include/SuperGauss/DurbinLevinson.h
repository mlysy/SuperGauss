/// @file DurbinLevinson.h

#ifndef DurbinLevinson_h
#define DurbinLevinson_h

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// FIXME: shouldn't happen whenever this file gets included
// using namespace Eigen; 

/// @brief Durbin-Levinson methods for Toeplitz matrices.
class DurbinLevinson {
 private:
  int N_; ///< Size of Toeplitz matrix.
  // temporary storage
  Eigen::VectorXd phi1_;
  Eigen::VectorXd phi2_;
  Eigen::VectorXd eps_;
  // the following vectors are dynamically reallocated only if their sizes change. 
  Eigen::VectorXd rx_;
  Eigen::VectorXd ry_;
  // Update the Durbin-Levinson coefficients.
  void update_phi(const Eigen::Ref<const Eigen::VectorXd>& acf, double nu, int ii);
 public:
  /// Constructor.
  DurbinLevinson(int N);
  /// Cross product with inverse Toeplitz.
  template <typename T1>
    double cross_prod(const Eigen::MatrixBase<T1>& M,
		      const Eigen::Ref <const Eigen::MatrixXd>& Xt,
		      const Eigen::Ref <const Eigen::MatrixXd>& Yt,
		      const Eigen::Ref <const Eigen::VectorXd>& acf, int calc_mode);
  /// Solve Toeplitz system.
  double solve(Eigen::Ref <Eigen::MatrixXd> X,
	       const Eigen::Ref <const Eigen::VectorXd>& acf,
	       const Eigen::Ref <const Eigen::MatrixXd>& Y);
  /// Multiply matrix by Cholesky of Toeplitz matrix or its inverse.
  void cholXZ(Eigen::Ref <Eigen::MatrixXd> Xt, Eigen::Ref<Eigen::MatrixXd> Zt,
	      const Eigen::Ref <const Eigen::VectorXd>& acf, bool ZtoX);
};

inline void DurbinLevinson::update_phi(const Eigen::Ref<const Eigen::VectorXd>& acf,
				       double nu, int ii) {
  phi2_.head(ii) = phi1_.head(ii).reverse();
  double rp = phi2_.head(ii).dot(acf.segment(1,ii));
  //rp = phi1_.head(ii).dot(acf.segment(1,ii).reverse());
  //phi2_.head(ii) = phi1_.head(ii).reverse();
  phi1_(ii) = (acf(ii+1) - rp)/nu;
  phi1_.head(ii) -= phi1_(ii) * phi2_.head(ii);
  return;
}


/// @param[in] N Size of Toeplitz matrix.
inline DurbinLevinson::DurbinLevinson(int N) {
  N_ = N;
  phi1_ = Eigen::VectorXd::Zero(N_);
  phi2_ = Eigen::VectorXd::Zero(N_);
  eps_ = Eigen::VectorXd::Zero(N_);
  rx_ = Eigen::VectorXd::Zero(1);
  ry_ = Eigen::VectorXd::Zero(1);
}

/// For `Tz = toeplitz(acf)`, calculates one of the following depending on the value of `calc_mode`:
///
/// - `calc_mode = 0`: `X' Tz^{-1} Y`
/// - `calc_mode = 1`: `X' Tz^{-1} X`
/// - `calc_mode = 2`: `trace(X' Tz^{-1} Y)`, where `X` and `Y` have the same dimensions.
///
/// @param[out] M Output matrix of size `d x k`, or vector of length `d = k` if `calc_mode = 2`.
/// @param[in] Xt Transpose of `X`: a matrix of size `d x N`.
/// @param[in] Yt Transpose of `Y`: a matrix of size `k x N`.  Ignored if `calc_mode = 1`.
/// @param[in] acf Vector of length `N` giving the first row/column of `Tz`.
/// @param[in] calc_mode Interger with values 0, 1, 2 determining the type of cross product to compute.
/// @return Log-determinant of `Tz` as a biproduct of the cross-product calculation.
template <typename T1>
inline double DurbinLevinson::cross_prod(const Eigen::MatrixBase<T1>& M,
					 const Eigen::Ref <const Eigen::MatrixXd>& Xt,
					 const Eigen::Ref <const Eigen::MatrixXd>& Yt,
					 const Eigen::Ref <const Eigen::VectorXd>& acf,
					 int calc_mode) {
  // required to use SelfAdjointView: https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
#define _M const_cast<Eigen::MatrixBase<T1>& >(M)
  int d,k;
  double nu, ldV;
  d = Xt.rows();
  k = (calc_mode != 1) ? Yt.rows() : 0;
  // resize rx_ and ry_ if necessary
  if(rx_.size() != d) rx_ = Eigen::VectorXd::Zero(d);
  if((calc_mode != 1) && (ry_.size() != k)) ry_ = Eigen::VectorXd::Zero(k);
  // initialize
  nu = acf(0);
  ldV = 0.0;
  _M.setZero();
  for(int ii = 0; ii < N_; ii++) {
    // variance
    if(ii > 0) nu *= (1.0-phi1_(ii-1)*phi1_(ii-1));
    ldV += log(nu);
    // residuals
    rx_ = Xt.leftCols(ii) * phi1_.head(ii).reverse();
    rx_ = Xt.col(ii) - rx_;
    if(calc_mode != 1) {
      ry_ = Yt.leftCols(ii) * phi1_.head(ii).reverse();
      ry_ = Yt.col(ii) - ry_;
    }
    // inner products
    if(calc_mode == 0) {
      // matrix inner product (default)
      rx_ /= nu;
      _M.noalias() += rx_ * ry_.transpose();
    } else if(calc_mode == 1) {
      // square norm (i.e. lower triangular only)
      rx_ /= sqrt(nu);
      _M.template selfadjointView<Eigen::Lower>().rankUpdate(rx_);
    } else {
      // elementwise inner product (i.e. diag(<X, Y>) )
      _M.array() += (rx_.array() * ry_.array())/nu;
    }
    // new coefficients
    if(ii < N_-1) {
      // coefficients
      update_phi(acf, nu, ii);
      // phi2_.head(ii) = phi1_.head(ii).reverse();
      // rp = phi2_.head(ii).dot(acf.segment(1,ii));
      // //rp = phi1_.head(ii).dot(acf.segment(1,ii).reverse());
      // //phi2_.head(ii) = phi1_.head(ii).reverse();
      // phi1_(ii) = (acf(ii+1) - rp)/nu;
      // phi1_.head(ii) -= phi1_(ii) * phi2_.head(ii);
    }
  }
  if(calc_mode == 1) {
    // upper triangular part of M
    _M.template triangularView<Eigen::Upper>() = M.transpose();
  }
  //std::cout << M << std::endl;
#undef _M
  return ldV;
}

/// Calculates `X = chol(Tz) Z` or `Z = chol(Tz)^{-1} X`, where `Tz = toeplitz(acf)`.
///
/// @param[in/out] Xt Transpose of `X`: a matrix of size `d x N`.
/// @param[in/out] Zt Transpose of `Z`: a matrix of size `d x N`.
/// @param[in] acf Vector of length `N` giving the first row/column of `Tz`.
/// @param[in] ZtoX Whether to multiply by the Cholesky factor (`ZtoX = false`) or its inverse (`ZtoX = true`).
void DurbinLevinson::cholXZ(Eigen::Ref <Eigen::MatrixXd> Xt, Eigen::Ref<Eigen::MatrixXd> Zt,
			    const Eigen::Ref <const Eigen::VectorXd>& acf, bool ZtoX) {
  double nu, sqrt_nu;
  nu = acf(0);
  // resize rx_ if necessary
  if(rx_.size() != Xt.rows()) rx_ = Eigen::VectorXd::Zero(Xt.rows());
  for(int ii = 0; ii < N_; ii++) {
    // variance
    if(ii > 0) nu *= (1-phi1_(ii-1)*phi1_(ii-1));
    sqrt_nu = sqrt(nu);
    // mean and actual conversion
    if(ii == 0) {
      rx_.setZero();
    } else {
      rx_ = Xt.leftCols(ii) * phi1_.head(ii).reverse();
    }
    if(ZtoX) {
      Xt.col(ii) = rx_ + sqrt_nu * Zt.col(ii);
    } else {
      Zt.col(ii) = (Xt.col(ii) - rx_)/sqrt_nu;
    }
    // new coefficients
    if(ii < N_-1) {
      // coefficients
      update_phi(acf, nu, ii);
      // coefficients
      // phi2_.head(ii) = phi1_.head(ii).reverse();
      // rp = phi2_.head(ii).dot(acf.segment(1,ii));
      // //rp = phi1_.head(ii).dot(acf.segment(1,ii).reverse());
      // //phi2_.head(ii) = phi1_.head(ii).reverse();
      // phi1_(ii) = (acf(ii+1) - rp)/nu;
      // phi1_.head(ii) -= phi1_(ii) * phi2_.head(ii);
    }
  }
  return;
}

/// Solves the system `Tz X = Y` for `Tz = toeplitz(acf)` using the Levinson-Trench-Zohar algorithm.
///
/// @param[out] Xt Transpose of `X`: a matrix of size `d x N`.
/// @param[in] acf Vector of length `N` giving the first row/column of `Tz`.
/// @param[in] Y Transpose of `Y`: a matrix of size `d x N`.
/// @return Log-determinant of `Tz` as a biproduct of the calculation.
inline double DurbinLevinson::solve(Eigen::Ref <Eigen::MatrixXd> Xt,
				    const Eigen::Ref <const Eigen::VectorXd>& acf,
				    const Eigen::Ref <const Eigen::MatrixXd>& Yt) {
  double xi, ldV;
  // resize rx_ = lambda if necessary
  if(rx_.size() != Xt.rows()) rx_ = Eigen::VectorXd::Zero(Xt.rows());
  // initialization
  phi1_(0) = 1.0; // phi1 = aa
  eps_(0) = acf(0);
  Xt.col(0) = Yt.col(0)/eps_(0);
  ldV = log(eps_(0));
  for(int ii=1; ii<N_; ii++) {
    // update xi and rx_
    phi2_.head(ii) = acf.segment(1,ii).reverse(); // needed for nocopy matmult
    xi = phi2_.head(ii).dot(phi1_.head(ii));
    rx_.noalias() = Xt.leftCols(ii) * phi2_.head(ii);
    xi /= -eps_(ii-1);
    rx_ = Yt.col(ii) - rx_;
    // Rprintf("xi = %f\n", xi);
    // std::cout << "lambda =\n" << rx_ << "\n" << std::endl;
    // update fwd coefficients
    phi2_.head(ii) = phi1_.head(ii).reverse(); // phi2 = bb = reverse phi1
    phi1_.segment(1, ii-1) += xi * phi2_.head(ii-1);
    phi1_(ii) = xi;
    // std::cout << "aa_ =\n" << phi1_.head(ii+1) << "\n" << std::endl;
    // update epsilon
    eps_(ii) = eps_(ii-1) * (1.0 - xi*xi);
    // update x
    rx_ /= eps_(ii);
    phi2_.head(ii) = phi1_.segment(1,ii).reverse();
    Xt.leftCols(ii).noalias() += rx_ * phi2_.head(ii).transpose();
    Xt.col(ii) = phi1_(0) * rx_;
    // update log-determinant
    ldV += log(eps_(ii));
  }
  return ldV;
}

#endif

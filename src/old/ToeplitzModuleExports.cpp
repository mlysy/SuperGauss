#include <Rcpp.h>
using namespace Rcpp;
#include "Toeplitz.h"

// wrapper class for Toeplitz objects
// this is necessary because the C++ Toeplitz class can't
// interact directly with R objects.
class Toeplitz_Rcpp {
private:
  Toeplitz *Tz_;
public:
  Toeplitz_Rcpp(int n) {
    Tz_ = new Toeplitz(n);
  }
  void setAcf(NumericVector acf) {
    Tz_->setAcf(REAL(acf));
    return;
  }
  NumericVector getAcf () {
    NumericVector acf(Tz_->size());
    Tz_->getAcf(REAL(acf));
    return acf;
  }
  NumericMatrix Multiply(NumericMatrix X) {
    int p = X.ncol();
    int n = X.nrow();
    NumericMatrix Y(n,p);
    for(int ii=0; ii<p; ii++) {
      Tz_->mult(&REAL(Y)[n*ii], &REAL(X)[n*ii]);
    }
    return Y;
  }
};

RCPP_MODULE(Toeplitz) {
  class_<Toeplitz_Rcpp>( "Toeplitz_cpp" )
    .constructor<double>()
    .method("setAcf", &Toeplitz_Rcpp::setAcf)
    .method("getAcf", &Toeplitz_Rcpp::getAcf)
    .method("Multiply", &Toeplitz_Rcpp::Multiply)
    ;
}

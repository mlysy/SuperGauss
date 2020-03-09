#include <Rcpp.h>
using namespace Rcpp;
#include "NormalToeplitz.h"

// R wapper functions to NormalToeplitz methods using XPtr

//[[Rcpp::export(".NormalToeplitz_constructor")]]
SEXP NormalToeplitz_constructor(int N) {
  NormalToeplitz *Nt = new NormalToeplitz(N);
  XPtr<NormalToeplitz> Nt_ptr(Nt, true);
  return Nt_ptr;
}

//[[Rcpp::export(".NormalToeplitz_logdens")]]
double NormalToeplitz_logdens(SEXP Nt_ptr, NumericVector z, NumericVector acf) {
  XPtr<NormalToeplitz> Nt(Nt_ptr);
  return Nt->logdens(REAL(z), REAL(acf));
}

//[[Rcpp::export(".NormalToeplitz_grad")]]
NumericVector NormalToeplitz_grad(SEXP Nt_ptr,
				  NumericVector z,
				  NumericMatrix dzdt, 
				  NumericVector acf,
				  NumericMatrix dadt,
				  int ntheta) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldt(ntheta);
	Nt->grad(REAL(dldt), REAL(z), REAL(dzdt), REAL(acf), REAL(dadt), ntheta);
	return dldt;
}

//[[Rcpp::export(".NormalToeplitz_hess")]]
NumericVector NormalToeplitz_hess(SEXP Nt_ptr,
				  NumericVector z,
				  NumericMatrix dzdt,
				  NumericMatrix d2zdt,
				  NumericVector acf,
				  NumericMatrix dadt,
				  NumericMatrix d2adt,
				  int ntheta) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericMatrix d2ldt(ntheta, ntheta);
	Nt->hess(REAL(d2ldt), REAL(z), REAL(dzdt), REAL(d2zdt), 
		 REAL(acf), REAL(dadt), REAL(d2adt), ntheta);
	return d2ldt;
}

//[[Rcpp::export(".NormalToeplitz_grad_full")]]
List NormalToeplitz_grad_full(SEXP Nt_ptr,
	NumericVector z, NumericVector acf) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	int N = Nt->size();
	NumericVector dldz(N);
	NumericVector dldacf(N);
	Nt->grad_full(REAL(dldz), REAL(dldacf), REAL(z), REAL(acf));
	List L = List::create(Named("dldz") = dldz, _["dldacf"] = dldacf);
	return L;
}

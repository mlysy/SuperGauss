#include <Rcpp.h>
using namespace Rcpp;
#include "NormalToeplitz.h"

// R wapper functions to NormalToeplitz methods using XPtr

//[[Rcpp::export(".NormalToeplitz_constructor")]]
SEXP NormalToeplitz_constructor(int n_, int p_) {
  NormalToeplitz *Nt = new NormalToeplitz(n_, p_);
  XPtr<NormalToeplitz> Nt_ptr(Nt, true);
  return Nt_ptr;
}

//[[Rcpp::export(".NormalToeplitz_logdens")]]
double NormalToeplitz_logdens(SEXP Nt_ptr, NumericVector z, NumericVector acf) {
  XPtr<NormalToeplitz> Nt(Nt_ptr);
  return Nt->logdens(REAL(z), REAL(acf));
}

//[[Rcpp::export(".NormalToeplitz_grad")]]
NumericVector NormalToeplitz_grad(SEXP Nt_ptr, NumericVector z, NumericVector dzdt, 
	NumericVector acf, NumericVector dacfdt) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldt(Nt->dim());
	Nt->grad(REAL(dldt), REAL(z), REAL(dzdt), REAL(acf), REAL(dacfdt));
	return dldt;
}

//[[Rcpp::export(".NormalToeplitz_hess")]]
NumericVector NormalToeplitz_hess(SEXP Nt_ptr,
	NumericVector z, NumericVector dzdt, NumericVector d2zdt,
	NumericVector acf, NumericVector dacfdt, NumericVector d2acfdt) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector d2ldt(Nt->dim(), Nt->dim());
	Nt->hess(REAL(d2ldt), REAL(z), REAL(dzdt), REAL(d2zdt), 
		REAL(acf), REAL(dacfdt), REAL(d2acfdt));
	return d2ldt;
}

//[[Rcpp::export(".NormalToeplitz_grad_full")]]
List NormalToeplitz_grad_full(SEXP Nt_ptr,
	NumericVector z, NumericVector acf) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldz(Nt->size());
	NumericVector dldacf(Nt->size());
	Nt->grad_full(REAL(dldz), REAL(dldacf), REAL(z), REAL(acf));
	List L = List::create(Named("dldz") = dldz, _["dldacf"] = dldacf);
	return L;
}

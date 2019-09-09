#include <Rcpp.h>
using namespace Rcpp;
#include "NormalToeplitz.h"

// R wapper functions to NormalToeplitz methods using XPtr

//[[Rcpp::export(".NormalToeplitz_constructor")]]
SEXP NormalToeplitz_constructor(int n, int p, bool hasToep) {
  NormalToeplitz *Nt = new NormalToeplitz(n, p, hasToep);
  XPtr<NormalToeplitz> Nt_ptr(Nt, true);
  return Nt_ptr;
}

//[[Rcpp::export(".NormalToeplitz_logdens")]]
double NormalToeplitz_logdens(SEXP Nt_ptr, NumericVector z, NumericVector acf) {
  XPtr<NormalToeplitz> Nt(Nt_ptr);
  return Nt->logdens(REAL(z), REAL(acf));
}

//[[Rcpp::export(".NormalToeplitz_logdens_")]]
double NormalToeplitz_logdens_(SEXP Nt_ptr, NumericVector z, Toeplitz* Tz) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	return Nt->logdens(REAL(z), Tz);
}

//[[Rcpp::export(".NormalToeplitz_grad")]]
NumericVector NormalToeplitz_grad(SEXP Nt_ptr,
	NumericVector z, NumericMatrix dzdt, 
	NumericVector acf, NumericMatrix dacfdt) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldt(Toep->p_);
	Nt->grad(REAL(dldt), REAL(z), REAL(dzdt), REAL(acf), REAL(dacfdt));
	return dldt;
}

//[[Rcpp::export(".NormalToeplitz_grad_")]]
NumericVector NormalToeplitz_grad_(SEXP Nt_ptr,
	NumericVector z, NumericMatrix* dzdt,
	Toeplitz* Tz, NumericMatrix dacfdt) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldt(Toep->p_);
	Nt->grad(REAL(dldt), REAL(z), REAL(dzdt), Tz, REAL(dacfdt));
	return dldt;
}

//[[Rcpp::export(".NormalToeplitz_hess")]]
NumericMatrix NormalToeplitz_hess(SEXP Nt_ptr,
	NumericVector z, NumericMatrix* dzdt, NumericVector d2zdt,
	NumericVector acf, NumericMatrix dacfdt, NumericVector d2acfdt) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericMatrix d2ldt(Toep->p_, Toep->p_);
	Nt->hess(REAL(d2ldt), REAL(z), REAL(dzdt), REAL(d2zdt), 
		REAL(acf), REAL(dacfdt), REAL(d2acfdt));
	return d2ldt;
}

//[[Rcpp::export(".NormalToeplitz_hess_")]]
NumericMatrix NormalToeplitz_hess_(SEXP Nt_ptr,
	NumericVector z, NumericMatrix* dzdt, NumericVector d2zdt,
	Toeplitz* Tz, NumericMatrix dacfdt, NumericVector d2acfdt) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericMatrix d2ldt(Toep->p_, Toep->p_);
	Nt->hess(REAL(d2ldt), REAL(z), REAL(dzdt), REAL(d2zdt),
		Tz, REAL(dacfdt), REAL(d2acfdt));
	return d2ldt;
}

//[[Rcpp::export(".NormalToeplitz_grad_full")]]
List NormalToeplitz_grad_full(SEXP Nt_ptr,
	NumericVector z, NumericVector acf) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldz(Toep->N_);
	NumericVector dldacf(Toep->N_);
	Nt->grad_full(REAL(dldz), REAL(dldacf), REAL(z), REAL(acf));
	List L = List::create(Named("dldz") = dldz, _["dldacf"] = dldacf);
	return L;
}

//[[Rcpp::export(".NormalToeplitz_grad_full_")]]
List NormalToeplitz_grad_full(SEXP Nt_ptr,
	Toeplitz* Tz, NumericVector acf) {
	XPtr<NormalToeplitz> Nt(Nt_ptr);
	NumericVector dldz(Toep->N_);
	NumericVector dldacf(Toep->N_);
	Nt->grad_full(REAL(dldz), REAL(dldacf), REAL(z), Tz);
	List L = List::create(Named("dldz") = dldz, _["dldacf"] = dldacf);
	return L;
}
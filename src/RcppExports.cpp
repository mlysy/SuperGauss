// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// DurbinLevinson_XZ
Eigen::MatrixXd DurbinLevinson_XZ(Eigen::MatrixXd X, Eigen::VectorXd acf);
RcppExport SEXP _SuperGauss_DurbinLevinson_XZ(SEXP XSEXP, SEXP acfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type acf(acfSEXP);
    rcpp_result_gen = Rcpp::wrap(DurbinLevinson_XZ(X, acf));
    return rcpp_result_gen;
END_RCPP
}
// DurbinLevinson_ZX
Eigen::MatrixXd DurbinLevinson_ZX(Eigen::MatrixXd Z, Eigen::VectorXd acf);
RcppExport SEXP _SuperGauss_DurbinLevinson_ZX(SEXP ZSEXP, SEXP acfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type acf(acfSEXP);
    rcpp_result_gen = Rcpp::wrap(DurbinLevinson_ZX(Z, acf));
    return rcpp_result_gen;
END_RCPP
}
// DurbinLevinson_Eigen
Rcpp::List DurbinLevinson_Eigen(Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::VectorXd acf, int calcMode);
RcppExport SEXP _SuperGauss_DurbinLevinson_Eigen(SEXP XSEXP, SEXP YSEXP, SEXP acfSEXP, SEXP calcModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type acf(acfSEXP);
    Rcpp::traits::input_parameter< int >::type calcMode(calcModeSEXP);
    rcpp_result_gen = Rcpp::wrap(DurbinLevinson_Eigen(X, Y, acf, calcMode));
    return rcpp_result_gen;
END_RCPP
}
// DurbinLevinson_Base
Rcpp::List DurbinLevinson_Base(NumericMatrix X, NumericMatrix Y, NumericVector acf, int calcMode);
RcppExport SEXP _SuperGauss_DurbinLevinson_Base(SEXP XSEXP, SEXP YSEXP, SEXP acfSEXP, SEXP calcModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    Rcpp::traits::input_parameter< int >::type calcMode(calcModeSEXP);
    rcpp_result_gen = Rcpp::wrap(DurbinLevinson_Base(X, Y, acf, calcMode));
    return rcpp_result_gen;
END_RCPP
}
// NormalToeplitz_constructor
SEXP NormalToeplitz_constructor(int N);
RcppExport SEXP _SuperGauss_NormalToeplitz_constructor(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(NormalToeplitz_constructor(N));
    return rcpp_result_gen;
END_RCPP
}
// NormalToeplitz_logdens
double NormalToeplitz_logdens(SEXP NTz_ptr, NumericVector z, NumericVector acf);
RcppExport SEXP _SuperGauss_NormalToeplitz_logdens(SEXP NTz_ptrSEXP, SEXP zSEXP, SEXP acfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type NTz_ptr(NTz_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    rcpp_result_gen = Rcpp::wrap(NormalToeplitz_logdens(NTz_ptr, z, acf));
    return rcpp_result_gen;
END_RCPP
}
// NormalToeplitz_grad
NumericVector NormalToeplitz_grad(SEXP NTz_ptr, NumericVector z, NumericMatrix dzdt, NumericVector acf, NumericMatrix dadt, int n_theta);
RcppExport SEXP _SuperGauss_NormalToeplitz_grad(SEXP NTz_ptrSEXP, SEXP zSEXP, SEXP dzdtSEXP, SEXP acfSEXP, SEXP dadtSEXP, SEXP n_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type NTz_ptr(NTz_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dzdt(dzdtSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dadt(dadtSEXP);
    Rcpp::traits::input_parameter< int >::type n_theta(n_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(NormalToeplitz_grad(NTz_ptr, z, dzdt, acf, dadt, n_theta));
    return rcpp_result_gen;
END_RCPP
}
// NormalToeplitz_hess
NumericMatrix NormalToeplitz_hess(SEXP NTz_ptr, NumericVector z, NumericMatrix dzdt, NumericMatrix d2zdt, NumericVector acf, NumericMatrix dadt, NumericMatrix d2adt, int n_theta);
RcppExport SEXP _SuperGauss_NormalToeplitz_hess(SEXP NTz_ptrSEXP, SEXP zSEXP, SEXP dzdtSEXP, SEXP d2zdtSEXP, SEXP acfSEXP, SEXP dadtSEXP, SEXP d2adtSEXP, SEXP n_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type NTz_ptr(NTz_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dzdt(dzdtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type d2zdt(d2zdtSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dadt(dadtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type d2adt(d2adtSEXP);
    Rcpp::traits::input_parameter< int >::type n_theta(n_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(NormalToeplitz_hess(NTz_ptr, z, dzdt, d2zdt, acf, dadt, d2adt, n_theta));
    return rcpp_result_gen;
END_RCPP
}
// NormalToeplitz_grad_full
List NormalToeplitz_grad_full(SEXP NTz_ptr, NumericVector z, NumericVector acf, bool calc_dldz, bool calc_dlda);
RcppExport SEXP _SuperGauss_NormalToeplitz_grad_full(SEXP NTz_ptrSEXP, SEXP zSEXP, SEXP acfSEXP, SEXP calc_dldzSEXP, SEXP calc_dldaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type NTz_ptr(NTz_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    Rcpp::traits::input_parameter< bool >::type calc_dldz(calc_dldzSEXP);
    Rcpp::traits::input_parameter< bool >::type calc_dlda(calc_dldaSEXP);
    rcpp_result_gen = Rcpp::wrap(NormalToeplitz_grad_full(NTz_ptr, z, acf, calc_dldz, calc_dlda));
    return rcpp_result_gen;
END_RCPP
}
// PCG_constructor
SEXP PCG_constructor(int n);
RcppExport SEXP _SuperGauss_PCG_constructor(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(PCG_constructor(n));
    return rcpp_result_gen;
END_RCPP
}
// PCG_Solve
NumericVector PCG_Solve(SEXP PCG_ptr, NumericVector acf, NumericVector y, double tol);
RcppExport SEXP _SuperGauss_PCG_Solve(SEXP PCG_ptrSEXP, SEXP acfSEXP, SEXP ySEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type PCG_ptr(PCG_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(PCG_Solve(PCG_ptr, acf, y, tol));
    return rcpp_result_gen;
END_RCPP
}
// Toeplitz_constructor
SEXP Toeplitz_constructor(int n);
RcppExport SEXP _SuperGauss_Toeplitz_constructor(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(Toeplitz_constructor(n));
    return rcpp_result_gen;
END_RCPP
}
// Toeplitz_setAcf
void Toeplitz_setAcf(SEXP Toep_ptr, NumericVector acf);
RcppExport SEXP _SuperGauss_Toeplitz_setAcf(SEXP Toep_ptrSEXP, SEXP acfSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Toep_ptr(Toep_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type acf(acfSEXP);
    Toeplitz_setAcf(Toep_ptr, acf);
    return R_NilValue;
END_RCPP
}
// Toeplitz_getAcf
NumericVector Toeplitz_getAcf(SEXP Toep_ptr);
RcppExport SEXP _SuperGauss_Toeplitz_getAcf(SEXP Toep_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Toep_ptr(Toep_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Toeplitz_getAcf(Toep_ptr));
    return rcpp_result_gen;
END_RCPP
}
// Toeplitz_Multiply
NumericMatrix Toeplitz_Multiply(SEXP Toep_ptr, NumericMatrix X);
RcppExport SEXP _SuperGauss_Toeplitz_Multiply(SEXP Toep_ptrSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Toep_ptr(Toep_ptrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Toeplitz_Multiply(Toep_ptr, X));
    return rcpp_result_gen;
END_RCPP
}
// Toeplitz_Solve
NumericMatrix Toeplitz_Solve(SEXP Toep_ptr, NumericMatrix X);
RcppExport SEXP _SuperGauss_Toeplitz_Solve(SEXP Toep_ptrSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Toep_ptr(Toep_ptrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Toeplitz_Solve(Toep_ptr, X));
    return rcpp_result_gen;
END_RCPP
}
// Toeplitz_Determinant
double Toeplitz_Determinant(SEXP Toep_ptr);
RcppExport SEXP _SuperGauss_Toeplitz_Determinant(SEXP Toep_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Toep_ptr(Toep_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Toeplitz_Determinant(Toep_ptr));
    return rcpp_result_gen;
END_RCPP
}
// Toeplitz_hasAcf
bool Toeplitz_hasAcf(SEXP Toep_ptr);
RcppExport SEXP _SuperGauss_Toeplitz_hasAcf(SEXP Toep_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Toep_ptr(Toep_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(Toeplitz_hasAcf(Toep_ptr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SuperGauss_DurbinLevinson_XZ", (DL_FUNC) &_SuperGauss_DurbinLevinson_XZ, 2},
    {"_SuperGauss_DurbinLevinson_ZX", (DL_FUNC) &_SuperGauss_DurbinLevinson_ZX, 2},
    {"_SuperGauss_DurbinLevinson_Eigen", (DL_FUNC) &_SuperGauss_DurbinLevinson_Eigen, 4},
    {"_SuperGauss_DurbinLevinson_Base", (DL_FUNC) &_SuperGauss_DurbinLevinson_Base, 4},
    {"_SuperGauss_NormalToeplitz_constructor", (DL_FUNC) &_SuperGauss_NormalToeplitz_constructor, 1},
    {"_SuperGauss_NormalToeplitz_logdens", (DL_FUNC) &_SuperGauss_NormalToeplitz_logdens, 3},
    {"_SuperGauss_NormalToeplitz_grad", (DL_FUNC) &_SuperGauss_NormalToeplitz_grad, 6},
    {"_SuperGauss_NormalToeplitz_hess", (DL_FUNC) &_SuperGauss_NormalToeplitz_hess, 8},
    {"_SuperGauss_NormalToeplitz_grad_full", (DL_FUNC) &_SuperGauss_NormalToeplitz_grad_full, 5},
    {"_SuperGauss_PCG_constructor", (DL_FUNC) &_SuperGauss_PCG_constructor, 1},
    {"_SuperGauss_PCG_Solve", (DL_FUNC) &_SuperGauss_PCG_Solve, 4},
    {"_SuperGauss_Toeplitz_constructor", (DL_FUNC) &_SuperGauss_Toeplitz_constructor, 1},
    {"_SuperGauss_Toeplitz_setAcf", (DL_FUNC) &_SuperGauss_Toeplitz_setAcf, 2},
    {"_SuperGauss_Toeplitz_getAcf", (DL_FUNC) &_SuperGauss_Toeplitz_getAcf, 1},
    {"_SuperGauss_Toeplitz_Multiply", (DL_FUNC) &_SuperGauss_Toeplitz_Multiply, 2},
    {"_SuperGauss_Toeplitz_Solve", (DL_FUNC) &_SuperGauss_Toeplitz_Solve, 2},
    {"_SuperGauss_Toeplitz_Determinant", (DL_FUNC) &_SuperGauss_Toeplitz_Determinant, 1},
    {"_SuperGauss_Toeplitz_hasAcf", (DL_FUNC) &_SuperGauss_Toeplitz_hasAcf, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_SuperGauss(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

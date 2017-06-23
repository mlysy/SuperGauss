// don't know how to do this properly for now...

#include <Rcpp.h>

void R_init_SuperGauss_modules(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}

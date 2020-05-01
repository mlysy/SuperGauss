/// @file mem_test.cpp
///
/// @brief `SuperGauss` tests which do not depend on R.
///
/// Used here to check for memory leaks with the LLVM [LeakSanitizer](https://clang.llvm.org/docs/LeakSanitizer.html).  To run:
///
/// ```
/// /usr/local/opt/llvm/bin/clang++ -I/Users/mlysy/Documents/R/SuperGauss/inst/include -I/usr/local/include -L/usr/local/lib -lfftw3 -fsanitize=address -g mem_test.cpp ; ASAN_OPTIONS=detect_leaks=1 ./a.out
/// ```

#include "VectorFFT.h"
#include "GSchur.h"
#include "Toeplitz.h"
#include "NormalToeplitz.h"
#include "PCG.h"

int main() {
  VectorFFT* fobj;
  fobj = new VectorFFT(100000);
  delete fobj;
  GSchur2K* gkobj;
  gkobj = new GSchur2K(100);
  delete gkobj;
  GSchurN* gnobj;
  gnobj = new GSchurN(100000);
  delete gnobj;
  Toeplitz* tobj;
  tobj = new Toeplitz(100000);
  delete tobj;
  NormalToeplitz* nobj;
  nobj = new NormalToeplitz(100000);
  delete nobj;
  PCG* pobj;
  pobj = new PCG(100000);
  delete pobj;
  return 0;
}

# SuperGauss 2.0.3

- Patch to support UCRT for R 4.2 on Windows.

# SuperGauss 2.0.0

## Major Changes

- The following functions are now **defunct**: `rSnorm()`, `dSnorm()`, `Snorm.grad()`, `Snorm.hess()`.  Replacements are `rnormtz()`, `rnormtz()` and the `$grad()` and `$hess()` methods in the `NormalToeplitz` class.

- Added `NormalToeplitz` and `NormalCirculant` classes in R and C++.

## Bug Fixes

- C++ library no longer calls `using namespace` in header files.

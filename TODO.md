## TODO

- [x] Fix R documentation for `Toeplitz` class.

- [x] Create `set_acf()` and `has_acf()` methods to `NormalToeplitz` and `NormalCirculant`, to make `logdens()`, `grad()`, etc are more efficient for multiple `z` with same `acf`.

- [x] Have `grad()`, etc. output the log-density as well.

- [x] Reimplement R-level `NormalToeplitz` convenience functions `dSnorm()`, `dSnorm.grad()`, etc.

	In fact, we'll only provide convenience functions `rnormtz()` and `dnormtz()`.  For the gradient and Hessian let's just use the `NormalToeplitz` class.
	
- [x] Add deprecation notes for `dSnorm()`, `rSnorm()`, etc.
	
- [x] Create class for Durbin-Levinson methods and add LTZ solve method.

- [x] Convert R `.` to `_` and rename other things.

- [x] Update vignette with new API.

- [x] Document `src/*Exports.cpp`.

- [x] Document PCG class.

- [x] Remove `using namespace` from header files.  Probably by wrapping the whole C++ library into a namespace.

	For simplicity just ended up using `Eigen::` to prefix its library members.

- [x] Finish `README.md`.

- [x] Update function documentation to use `man-roxygen` and `examples`.  For the former, it means using the `@template` tags extensively for **roxygen2** documentation (don't forget to add `man-roxygen` to `.Rbuildignore`).  For the latter, it means that examples for R function `foo` should be in plain R format (i.e., no `'#`) in `examples/foo.R`, and referenced in the documentation via `@example examples/foo.R`.  The point is that it's much easier to write examples directly in R than have to comment them out all the time.  Also, you can share examples between functions this way, etc.

- [x] Need to `.gitignore` some files and/or hide things in `tests/dontrun`.

- [x] ~~**Set up as R package** *(Yun)*~~.  `SuperGauss` should be the name of the **R** package.  Merge the Durbin-Levinson functions back in as they are required for some variants of `rSnorm`, and for Cholesky decomposition in general.  Use `.Rbuildignore` and `.gitignore` to ignore things that shouldn't be in the **R** package (such as `setup.R`) and the Git repo (such as `.o`, `.dll`, `.Rhistory`), respectively.  Note that **R** by default doesn't include anything in folders named `old` or any file starting with ".". **Edit:** `Toeplitz` folder should be removed, such that top-level of `SuperGaus` is already the package.

- [x] ~~**Roxygen documentation** *(Martin)*~~.  Write the **R** function documentation and make sure it is `devtools`-compatible.

    In documentation, `@title` is for very short description of function.  `@description` gives a bit more detail.  Extra detail, (e.g., acf definition, restrictions on inputs) are given in `@details`.  Be *consistent*: for example:
	
    - don't sometimes use ACF, sometimes autocorrelation (should probably use the longer one all the time, same for MSD.  Otherwise, write it once per function then (ACF) right after).  	
    - Consistency of capital letters and punctuation in titles.  Let's say for all documentation: Capitalize first word only, period at the end.
    - Consistency of `code` vs latex notation.  Typically, avoid latex, unless `code` becomes too confusing.  So most of the time I only use latex in `\deqn`.
    - Let's define all `acf` inputs as `acf` (this goes against general rule of not naming function arguments by built-in functions).


- [x] ~~**Toeplitz class documentation** *(Yun & Martin)*~~.  Set up `Toeplitz-Class.R` to do this *(Martin)*.  Then documentation is exactly like regular functions *(Yun)*.

- [x] **Package vignette** *(Yun)*.  This is basically a quick tutorial for the package and what it can do.  Suggestion: fBM + stochastic drift.  Simulate data, use Newton-Raphson to recover parameters, simulate posterior dstr of drift using Dietrich algorithm.  The first section should be instructions on obtaining `fftw`.  Write this using `Rmarkdown` as per the following instructions: <http://r-pkgs.had.co.nz/vignettes.html>

- [x] ~~**Header-only library** *(Yun & Martin)*~~.  Convert all **C++** code to `.h` files *(Yun)*.  Put these in `inst/include/` and make sure the **R** wrappers know where to look *(Martin)*. 

- [x] ~~**More formal testing** *(Yun)*~~.  Set up some tests for all the functions using the `testthat` package as per the following instructions: <http://http://r-pkgs.had.co.nz/tests.html>.  This completely automates the test procedure, such that anytime something in the package or **C++** code gets modified, or if we want to see if things work on linux or Mac, we can essentially do this with a single line of code.

- [x] **Modify list of provided acf's** *(Yun)*.  There are too many acf's right now, and now enough explanation for the unfamiliar ones.  Let's reduce this to the following: `fbm.acf`, `pex.acf`, i.e., power exponential: `pex.acf(N, dt, lambda, rho) = exp(-abs((1:N-1)*dt/lambda)^rho)`.  Let's use the following functions:

    ```r
    fbm.msd(tseq, alpha)
	pex.acf(tseq, lambda, rho)
	matern.acf(tseq, lambda, nu)
	msd2acf(msd) # only gives correct result if msd is equally spaced
	acf2msd(acf)
	acf2incr(acf)
    ```

    Then in the vignette we can define

    ```r
    fbm.acf <- function(alpha, N, dt) {
      msd2acf(msd = fbm.msd((1:N-1)*dt, alpha))
    }
    ```

- [x] **Rename `Toep.mult`** *(Yun)*.  Let's do: `toep.mult <- function(acf, x)`

- [x] **Finalize vignette** *(Martin + Yun)*.  Add a bit more details to vignette, e.g. what is the `fbm.acf`?  Also nice to add CI for one of the calculations.  Appendix: installation notes.

- [x] **Optim issue** *(Martin)*. In the vignette, must document the following issue.  Probably in other places in the package too

    ```r
    N <-100
    myacf <- function(lambda) exp(-(1:N-1)/lambda)
    Toep <- Toeplitz(n = N)

    loglik <- function(lambda) {
      acf <- Toep$setAcf(myacf(lambda))
      dSnorm(X, acf = acf, log = TRUE)
    }

    fit <- optim(fun = loglik, par = 1, control = list(fnscale = -1))
    # original value of Toep is now overwritten!
    ```

- [x] **Optimize code for `Snorm.grad` and `Snorm.hess`** *(Martin)*.  Perhaps fewer for-loops, fewer `apply`s, etc.

- [x] Separate C++ `Toeplitz` class from `GSchur` class, and `NormalToeplitz` class.  The point is that `Toeplitz` is templated on either `GSchur` or `LTZ`.

    As for `NormalToeplitz`, here's a potential class design:
	
	```c
	class NormalToeplitz {
	  private:
	  Toeplitz *Tz_; ///> Internal pointer to Toeplitz matrix.
	                 ///> Always postfix internal variables with "_"
	  public:
	  /// Constructor.
	  NormalToeplitz(int N);
	  /// Perhaps also need a ctor which doesn't allocate the Toeplitz memory internally?
	  /// Make sure that all methods take `N` is taken from `Tz_`, 
	  /// i.e., don't use private `N_` because external memory version won't know what this is.
	  NormalToeplitz(void);
	  /// Log-Density.
	  /// TODO: should we use std::vector<double> instead? 
	  /// is there a speed difference?
	  double logdens(const double* z, const double* acf) {
	    Tz_->setAcf(acf);
	    return logdens(z, Tz_);
	  }
	  /// this is the preferred form for R, 
	  /// because we don't want to reallocate memory for every call.
	  /// note that `Tz` can't be `const` because we potentially modify its internal structure.
	  double logdens(const double* z, Toeplitz* Tz);
	  /// Full gradient.
	  /// The outputs `dldz` and `dldacf` are each the length of `z`.
	  void grad(double* dldz, double* dldacf,
	            const double* z, const double* acf);
	  /// Gradient with respect to theta.
	  /// The output `dldt` is the same length as `theta`.
	  void grad(double* dldt,
	           const double* z, const double* dzdt, 
	           const double* acf, const double* dacfdt);
	  /// Same thing with Toeplitz input
	  void grad(double* dldt,
	           const double* z, const double* dzdt, 
	           Toeplitz* Tz, const double* dacfdt);
	  /// Hessian with respect to theta.
	  /// Makes no sense to have full hessian since it scales as `O(N^2)`.
	  /// Names!
	  void hess(double* d2ldt,
	            const double* z, const double* dzdt, const double* d2zdt, 
	            const double* acf, const double* dacfdt, const double* d2acfdt);
	  void hess(double* d2ldt,
	            const double* z, const double* dzdt, const double* d2zdt, 
	            Toeplitz* Tz, const double* dacfdt, const double* d2acfdt);
	};
	```


- [x] Implement gradient in Stan.  We will need to talk to the Stan development team about this.  Here are some of the links I found for help:

    - [Full reference](https://arxiv.org/abs/1509.07164)
	- [Simple example]( https://github.com/stan-dev/math/wiki/Adding-a-new-function-with-known-gradients).  What this is missing however is memory allocation, i.e., you do not want to allocate memory every time Stan computes the gradient!
	
	My suggestion is to write to the stan developers via the [Stan forum](https://discourse.mc-stan.org/).  We can tell them that our C++ code works something like this:
	
	```c
	#include <NormalToepliz.h>
	
	// Create an object once that does all memory allocation
	NormalToeplitz nt(N); // int N = size of problem
	
	// Compute loglikelihood
	double ll = nt.loglik(z, acf);
	// Compute loglikelihood and gradient wrt z and acf simultaneously
	nt.lgrad(ll, dldz, dldacf);
	```

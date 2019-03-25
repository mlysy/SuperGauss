## TODO

- [ ] Update function documentation to use `man-roxygen` and `examples`.  For the former, it means using the `@template` tags extensively for **roxygen2** documentation (don't forget to add `man-roxygen` to `.Rbuildignore`).  For the latter, it means that examples for R function `foo` should be in plain R format (i.e., no `'#`) in `examples/foo.R`, and referenced in the documentation via `@example examples/foo.R`.  The point is that it's much easier to write examples directly in R than have to comment them out all the time.  Also, you can share examples between functions this way, etc.

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

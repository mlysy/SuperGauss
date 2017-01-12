### TODO

1. **Package vignette** *(Yun)*.  This is basically a quick tutorial for the package and what it can do.  Suggestion: fBM + stochastic drift.  Simulate data, use Newton-Raphson to recover parameters, simulate posterior dstr of drift using Dietrich algorithm.  The first section should be instructions on obtaining `fftw`.  Write this using `Rmarkdown` as per the following instructions: <http://r-pkgs.had.co.nz/vignettes.html>
2. **Roxygen documentation** *(Martin)*.  Write the **R** function documentation and make sure it is `devtools`-compatible.
3. **Header-only library** *(Yun & Martin)*.  Convert all **C++** code to `.h` files *(Yun)*.  Put these in `inst/include/` and make sure the **R** wrappers know where to look *(Martin)*. 

# SuperGauss: Superfast Likelihood Inference for Stationary Gaussian Time Series

*Yun Ling, Martin Lysy*

---

### Description

Likelihood evaluations for stationary Gaussian time series are typically obtained via the Durbin-Levinson algorithm, which scales as O(N^2) in the number of time series observations.  This package provides a "superfast" O(N log^2 N) algorithm written in C++, crossing over with Durbin-Levinson around N = 300.  Efficient implementations of the score and Hessian functions are also provided, leading to superfast versions of inference algorithms such as Newton-Raphson and Hamiltonian Monte Carlo.  The C++ code provides a `Toeplitz` matrix class and a `NormalToeplitz` distribution class packaged as a header-only library, to simplify low-level usage in other packages and outside of R.  A complete description of the algorithms is available in this [preprint](doc/SuperGauss_preprint.pdf).

### Installation

Version 1.0.2 is available on [CRAN](https://CRAN.R-project.org/package=SuperGauss), but in order to install this much-updated development version from GitHub, it is necessary to first install the [FFTW](http://www.fftw.org/) library:

- **Windows**:

	- Download the precompiled FFTW binary from [here](http://www.fftw.org/install/windows.html) for your version of Windows (32 or 64 bit).
   - In order to install **SuperGauss** directly as per the instructions below, you **must** unzip the DLLs to a folder called `C:/fftw`, and add this location to the system `PATH` variable, as explained [here](https://www.java.com/en/download/help/path.xml).
   - If for some reason you can't install to `C:/fftw`, then in the **SuperGauss** source folder, open the file `src/Makevars.win` and replace the instances of `C:/fftw` with the location where the **FFTW** library is installed.  Make sure you modify the system `PATH` variable accordingly.

- **macOS**:

	- Download the latest version of FFTW from [here](http://www.fftw.org/download.html).
   - Install the library yourself by opening Terminal in the downloaded folder, then:

		```bash
		./configure
		make
		make install
		```

		Further instructions to customize the installation are available [here](http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix).
		
- **Linux:** Should be straightforward.

- **Testing:** To make sure the FFTW library is correctly installed, try building the R package  [**fftw**](https://CRAN.R-project.org/package=fftw) from source.  

Once you have successfully installed FFTW, **SuperGauss** can be installed from GitHub via the command:

```r
devtools::install_github("mlysy/SuperGauss", ref = "devel_martin")
```

### C++ API

The `NormalToeplitz` class corresponds to the distribution 

```
z ~ N( 0, Tz = Toeplitz(acf) ),
```

where `z` is an `N`-dimensional vector and `Tz` is a symmetric positive-definite Toeplitz variance matrix with "autocorrelation" (i.e., first row/column) `acf`.  Extensive member documentation is provided by Doxygen code comments, but very briefly we have the following:

```c
NormalToeplitz NTz(N); // constructor

// log-density.  both z and acf are `double*`
ld = NTz.logdens(z, acf); 

// log-density gradient wrt z and/or acf
// dldz, dlda: `double*` outputs
// z, acf: `double*` inputs
// calc_dldz, calc_dlda: `bool` to calculate one or both gradients.
// when e.g., `calc_dldz = false`, `dldz` is not modified.
NTz.grad_full(dldz, dlda, z, acf, calc_dldz, calc_dlda)
```

### R API

TBD

### TODO

- [x] Fix R documentation for `Toeplitz` class.

- [x] Create `set_acf()` and `has_acf()` methods to `NormalToeplitz` and `NormalCirculant`, to make `logdens()`, `grad()`, etc are more efficient for multiple `z` with same `acf`.

- [x] Have `grad()`, etc. output the log-density as well.

- [x] Reimplement R-level `NormalToeplitz` convenience functions `dSnorm()`, `dSnorm.grad()`, etc.

	In fact, we'll only provide convenience functions `rnormtz()` and `dnormtz()`.  For the gradient and Hessian let's just use the `NormalToeplitz` class.
	
- [x] Create class for Durbin-Levinson methods and add LTZ solve method.

- [x] Convert R `.` to `_` and rename other things.

- [x] Update vignette with new API.

- [ ] Document `src/*Exports.cpp`.

- [ ] Document PCG class.

- [ ] Remove `using namespace` from header files.  Probably by wrapping the whole C++ library into a namespace.

- [ ] Finish this `README.md`.

#--- setup for Toeplitz package -------------------------------------------------

# for recompiling package
# first quit R, then setwd() to where setup.R is found. then:
# pkg.path <- "D:/GitHub/SuperGauss"
pkg.path <- "c:/Users/Jerome/Documents/R/SuperGauss"

require(Rcpp)
require(devtools)

# regenerates Rcpp interface (i.e., RcppExports)
compileAttributes(pkgdir = pkg.path)
document(pkg = pkg.path)
install(pkg = pkg.path) # installs the package
build(pkg = pkg.path) # builds a tar.gz file

# restart R before testing changes
testthat::test_package("SuperGauss")

build(pkg = pkg.path)

# build windows binary
# NOTE: this will also install the package
build(binary = TRUE)

pkg.path <- "c:/Users/Jerome/Documents/R/test/fftw"
cmd <- file.path(R.home(component = "bin"),
                 paste0("R CMD INSTALL ", pkg.path))
compiled <- system(cmd)

# First Time Installation -------------------------------------------------
# setwd("D:/GitHub/SuperGauss")
require(Rcpp)
require(RcppEigen)

pkg.name <- "Toeplitz"
pkg.path <- getwd()

# for safest results, delete ALL traces of previous package versions
remove.packages(pkg.name)
unlink(file.path(pkg.path, pkg.name), recursive = TRUE)

# minimum number of files to get the package started
Rcpp.package.skeleton(name = pkg.name, path = pkg.path, cpp_files =
                        c("Toeplitz.h", "Toeplitz.cpp", "VectorFFT.cpp", "VectorFFT.h", "GSchur.cpp", "GSchur.h", "Wrapper.cpp"),
                      example_code = FALSE, force = TRUE, module = TRUE)

## removing the generated Sample code
rm <- file.remove(file.path(pkg.path, pkg.name, "src", "Num.cpp"))
rm <- file.remove(file.path(pkg.path, pkg.name, "src", "rcpp_module.cpp"))
rm <- file.remove(file.path(pkg.path, pkg.name, "src", "stdVector.cpp"))

cat("loadModule(\"Toeplitz\", TRUE)\n", file = file.path(pkg.path, pkg.name, "R", "zzz.R"))

# link to RcppEigen
DESCRIPTION <- read.dcf(file = file.path(pkg.path, pkg.name, "DESCRIPTION"))
DESCRIPTION[,"LinkingTo"] <- paste0(DESCRIPTION[,"LinkingTo"], ", RcppEigen")
write.dcf(DESCRIPTION, file = file.path(pkg.path, pkg.name, "DESCRIPTION"))
compileAttributes(pkgdir = file.path(pkg.path, pkg.name))

# link to FFTW library (windows only)
cat("PKG_CXXFLAGS=-I\"C:/fftw\" `${R_HOME}/bin/Rscript -e \"Rcpp:::CxxFlags()\"`",
    "PKG_LIBS=-L\"C:/fftw\" -lfftw3-3 `${R_HOME}/bin/Rscript -e \"Rcpp:::LdFlags()\"`" ,
    sep = "\n", file = file.path(pkg.path, pkg.name, "src", "Makevars.win"))

# install package
cmd <- file.path(R.home(component = "bin"), paste0("R CMD INSTALL ", pkg.name))
compiled <- system(cmd)


#  ------------------------------------------------------------------------

setwd("D:/GitHub/SuperGauss/tests/testthat")
require(testthat)
require(SuperGauss)
source("SuperGauss-test-functions.R")

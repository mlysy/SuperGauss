#--- setup for Toeplitz package -------------------------------------------------

# for recompiling package
# quit R, then setwd() to where setup.R is found. then:
# setwd("D:/GitHub/SuperGauss")

require(Rcpp)
require(devtools)

compileAttributes() # regenerates Rcpp interface (i.e., RcppExports)
document()
install() # installs the package
# build() # builds a tar.gz file

# restart R before testing changes


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

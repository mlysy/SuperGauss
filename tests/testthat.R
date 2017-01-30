library(testthat)
library(SuperGauss)
library(RUnit)

test_check("SuperGauss")

# dimension
N <- round(abs(rnorm(n = 1, mean = 200, sd = 10)))
d <- round(abs(rnorm(n = 1, mean = 10, sd = 3)))

# creating the class
Toep <- new(Toeplitz, N)


# dimension check for Toeplitz
checkEquals(Toep$DimCheck(), N)


# input acf
lambda <- 100
acf <- exp(-(1:N)^2/lambda) * lambda

# input data can be vector, matrix or array, only requirement is total length should be N
Toep$AcfInput(acf)
Toep$AcfInput(matrix(acf, N, 1))
Toep$AcfInput(array(acf, c(N, 1, 1)))

Toep.mat <- toeplitz(acf)

# if total length is wrong, it will stop, but the previous information is not overwritten
Toep$AcfInput(c(acf, 0))
Toep$AcfInput(acf[-1])

# Toeplitz - matrix multiplication
X <- matrix(rnorm(N * d), N, d)

checkEquals(Toep$Mult(X), Toep.mat %*% X)

# inverse Toeplitz - matrix multiplication

checkEquals(Toep.mat %*% Toep$Solve(X), X)

# determinant 

checkEquals(log(det(Toep.mat)), Toep$Det())

# traceProd


# traceDeriv


# dSnorm


# Snorm.grad


# Snorm.hess

## some idea:
## we should check the accuracy of GSchur againts Durbin-Levinson or Choleski?
## if we wannna compare the outcome of 
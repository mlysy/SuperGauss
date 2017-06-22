require(SuperGauss)

N <- 20
Toep <- Toeplitz(N)

lambda <- .1
dT <- 1/60
acf <- exp2.acf(lambda = lambda, dT = dT, N = N, incr = TRUE)
Toep$setAcf(acf)

#--- det test ------------------------------------------------------------------

determinant(Toep)
2*sum(log(diag(chol(toeplitz(acf)))))
SuperGauss:::DurbinLevinson_Eigen(X = matrix(1, N), Y = matrix(1, N), acf = acf, calcMode = 1)$ldV

#--- solve test ----------------------------------------------------------------

solveV <- function(V, x) {
  C <- chol(V)
  backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
}


x <- matrix(rnorm(N))
y <- toeplitz(acf) %*% x

C <- chol(toeplitz(acf))

# test
n <- 10
range(solve(toeplitz(acf[1:n]), y[1:n]) - solveV(V = toeplitz(acf[1:n]), x = y[1:n]))

range((solveV(toeplitz(acf), y) - x)/abs(x))
range((solve(Toep, y) - x)/abs(x))

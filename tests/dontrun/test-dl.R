require(SuperGauss)

n <- 10
p <- 2
acf <- exp(-(1:n))
X <- matrix(rnorm(n*p), n, p)
Y <- matrix(rnorm(n*p), n, p)

SuperGauss:::DurbinLevinson_Eigen(X = X, Y = Y, acf = acf, calcMode = 0)$IP
crossprod(X, solve(toeplitz(acf), Y))


SuperGauss:::DurbinLevinson_Eigen(X = X, Y = matrix(0),
                                  acf = acf, calcMode = 1)$IP
crossprod(X, solve(toeplitz(acf), X))

SuperGauss:::DurbinLevinson_Eigen(X = X, Y = Y, acf = acf, calcMode = 2)$IP
diag(crossprod(X, solve(toeplitz(acf), Y)))

SuperGauss:::DurbinLevinson_Eigen(X = X, Y = Y, acf = acf, calcMode = 2)$ldV
determinant(toeplitz(acf))$mod

N <- 10
p <- 2
acf <- exp(-(1:N - 1))

Z <- matrix(rnorm(N * p), N, p)
cholZX(Z = Z, acf = acf) - (t(chol(toeplitz(acf))) %*%  Z)

X <- matrix(rnorm(N * p), N, p)
cholXZ(X = X, acf = acf) - solve(t(chol(toeplitz(acf))), X)

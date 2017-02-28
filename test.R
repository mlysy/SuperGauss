#--- basic test of SuperGauss package --------------------------------------

require(SuperGauss)

# Toeplitz class
tr <- function(X) sum(diag(X)) # matrix trace
n <- 200
d <- 4
acf <- exp(-(1:n)^2/102)*102
x <- matrix(rnorm(n*d), n, d)
y <- matrix(rnorm(n*d), n, d)
z <- matrix(rnorm(n*d), n, d)
T1 <- new(Toeplitz, n)
T1$AcfInput(acf)
## dimension check
T1$DimCheck()
## InverseProd part
msg <- T1$Mult(x)
range(msg - toeplitz(acf) %*% x)
msg <- T1$Solve(y)
msg1 <- T1$Solve(x)
msg2 <- T1$Solve(z)
range(toeplitz(acf) %*% msg - y)
# determinant part
T1$Det() - log(det(toeplitz(acf)))
# determinant part
acf2 <- exp(-(1:n)^2/12)*12
T1$AcfInput(acf2)
T1$Det() - log(det(toeplitz(acf2)))

#--- traceProd test --------------------------------------------------------

require(SuperGauss)
tr <- function(X) sum(diag(X)) # matrix trace
n <- 46
acf <- exp(-(1:n)^2/102)*102
T1 <- new(Toeplitz, n)
T1$AcfInput(acf)
acf2 <- rnorm(n)
acf3 <- rnorm(n)
acf4 <- rep(0, n)
T1$TraceProd(acf2) - tr(solve(toeplitz(acf)) %*% toeplitz(acf2))
T1$TraceProd(acf3) - tr(solve(toeplitz(acf)) %*% toeplitz(acf3))
T1$TraceProd(acf4) - tr(solve(toeplitz(acf)) %*% toeplitz(acf4))
## delete
rm(T1)


#--- traceDerv test --------------------------------------------------------

require(SuperGauss)
require(numDeriv)
tr <- function(X) sum(diag(X)) # matrix trace
# derivative trace in R
traceDerv.R <- function(acf, acf2, acf3, debug = F){
  if(debug) browser()
  T1 <- toeplitz(acf)
  T2 <- toeplitz(acf2)
  T3 <- toeplitz(acf3)
  n <- length(acf)
  sT1 <- solve(T1)
  phi <- matrix(sT1[,1], n, 1)
  phi2 <- - sT1 %*% T2 %*% phi

  phi <- as.vector(phi)
  phi2 <- as.vector(phi2)

  L1.phi2 <- toeplitz(phi2)
  L1.phi2[upper.tri(L1.phi2)] <- 0

  L1.phi <- toeplitz(phi)
  L1.phi[upper.tri(L1.phi)] <- 0

  L2.phi2 <- toeplitz(c(0, rev(phi2[-1])))
  L2.phi2[upper.tri(L2.phi2)] <- 0

  L2.phi <- toeplitz(c(0, rev(phi[-1])))
  L2.phi[upper.tri(L2.phi)] <- 0

  trace <- - phi2[1] * tr(sT1 %*% T3)

  trace <- trace + 2 * tr(L1.phi2 %*% t(L1.phi) %*% T3)

  # knowing that T3 = 1/acf3[1] * (L3 L3' - L4 L4')
  L3 <- toeplitz(acf3)
  L3[upper.tri(L3)] <- 0
  L4 <- toeplitz(c(0, acf3[-1]))
  L4[upper.tri(L4)] <- 0

  trace <- trace - 2 * tr(L2.phi2 %*% t(L2.phi) %*% T3)

  -trace / phi[1]
}

# acf and acf1 cannot be arbitary vectors
# test, acf = exp(-(i/theta)^lambda)
acf.fun0 <- function(theta, lambda, N){
  acf <- 1:N
  acf <- -(acf / theta)^lambda
  exp(acf)
}

# d.acf = lambda * exp(-(i/theta)^lambda) * (i/theta)^lambda / theta
acf.fun1 <- function(theta, lambda, N){
  acf <- 1:N
  acf <- (acf / theta)^lambda
  acf <- lambda * exp(-acf) * acf / theta
  acf
}

n <- 6
acf <- exp(-(1:n)^2/102)*102
T1 <- new(Toeplitz, n)
T1$AcfInput(acf)
acf2 <- rnorm(n)
acf3 <- rnorm(n)
T1$TraceDeriv(acf2, acf3) - traceDerv.R(acf, acf2, acf3)

range(grad(func = acf.fun0, x = 2, lambda = 2, N = 1) - acf.fun1(2, 2, 1))
# good

N <- 46
theta <- 2
lambda <- 1
acf <- acf.fun0(theta, lambda, N)
acf1 <- acf.fun1(theta, lambda, N)
acf2  <- rnorm(N)

T1 <- toeplitz(acf)
sT1 <- solve(T1)
T2 <- toeplitz(acf1)
T3 <- toeplitz(acf2)

traceDerv.R(acf, acf1, acf2) - tr(sT1 %*% T2 %*% sT1 %*% T3)
# good
## delete
gc()

#--- Durbin-Levinson functions ---------------------------------------------

# inner product
n <- 50
d <- 5
k <- 3
acf <- exp(-(1:n)^2/102)*102
T <- toeplitz(acf)

X <- matrix(rnorm(n*d), n, d)
Y <- matrix(rnorm(n*k), n, k)

# calcMode == 1
IP <- crossprod(X, solve(T, X))
ldV <- determinant(T, log = TRUE)$mod[1]
tmp <- DurbinLevinsonEigen(X = X, Y = as.matrix(0),
                           acf = acf, calcMode = 1)
range(tmp$IP - IP)
abs(tmp$ldV - ldV)

# calcMode == 0
IP <- crossprod(X, solve(T, Y))
ldV <- determinant(T, log = TRUE)$mod[1]
tmp <- DurbinLevinsonEigen(X = X, Y = Y,
                           acf = acf, calcMode = 0)
range(tmp$IP - IP)
abs(tmp$ldV - ldV)

# calcMode == 2
IP <- diag(crossprod(X, solve(T, Y)))
ldV <- determinant(T, log = TRUE)$mod[1]
tmp <- DurbinLevinsonEigen(X = X, Y = Y,
                           acf = acf, calcMode = 2)
range(tmp$IP - IP)
abs(tmp$ldV - ldV)

# X -> Z

n <- 50
d <- 5
acf <- exp(-(1:n)^2/102)*102
T <- toeplitz(acf)

# X -> Z
X <- matrix(rnorm(n*d), n, d)
Z <- solve(t(chol(T)), X)
Z2 <- toeplitzXZ(X = X, acf = acf)
range(Z-Z2)

# Z -> X
Z <- matrix(rnorm(n*d), n, d)
X <- t(chol(T)) %*% Z
X2 <- toeplitzZX(Z = Z, acf = acf)
range(X-X2)

# time check
nrep <- 100
n <- 2000
k <- 100
# 1
Toep1 <- new(Toeplitz, n)
system.time({
  for(ii in 1:nrep){
    acf <- exp(-(1:n)^2/102)*102
    X <- matrix(rnorm(n*k), n, k)
    Toep1$AcfInput(acf)
    Toep1$Mult(X)
  }
})

#--- depreciated -----------------------------------------------------------

if(FALSE) {
# 2
Toep2 <- new(Toeplitz, n, k)
system.time({
  for(ii in 1:nrep){
    acf <- exp(-(1:n)^2/102)*102
    X <- matrix(rnorm(n*k), n, k)
    Toep2$AcfInput(acf)
    Toep2$Mult(X)
  }
})
}

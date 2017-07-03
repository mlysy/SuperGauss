# decomposition of Toeplitz matrix with 0 on diagonal

# generate non-symmetric toeplitz matrix.  xr and xc are the first row and column, and xr[1] == xc[1].
toeplitz2 <- function(xr, xc) {
  if((length(xr) != length(xc)) || (xr[1] != xc[1])) stop("Incorrect row and column specification.")
  n <- length(xr)
  ans <- toeplitz(xr)
  ans[lower.tri(ans)] <- toeplitz(xc)[lower.tri(ans)]
  ans
}

N <- 5
gam <- rnorm(N)
U1 <- toeplitz2(gam, c(gam[1], rep(0, N-1)))
U2 <- toeplitz2(c(0, gam[-1]), rep(0, N))
(tcrossprod(U1) - tcrossprod(U2))/gam[1] - toeplitz(gam)

N <- 5
gam <- c(0, rnorm(N-1))
U1 <- toeplitz2(gam, c(gam[1], rep(0, N-1)))
U2 <- toeplitz2(c(0, gam[-1]), rep(0, N))
(tcrossprod(U1) - tcrossprod(U2))

#--------------------------------------------------------------------------

require(SuperGauss)

N <- 6
acf <- exp(-(1:N)^2*runif(1))
Toep <- Toeplitz(acf = acf)
iT <- solve(Toep)
phi <- SuperGauss:::.Toeplitz_getPhi(Toep$cpp_ptr)

# ok let's try GS decomp
sum(diag(iT))
L1 <- toeplitz2(c(phi[1], rep(0, N-1)), phi)
L2 <- toeplitz2(rep(0, N), c(0, phi[N:2]))

(tcrossprod(L1) - tcrossprod(L2))/phi[1]
iT

sum((N:1) * (phi^2 - c(0, phi[N:2])^2))/phi[1]

N*phi[1] + sum((N-1):1 * (phi[2:N]^2 - phi[N:2]^2)/phi[1])
sum(seq(N, by = -2, len = N) * phi[1:N]^2)/phi[1]
sum(diag(iT))

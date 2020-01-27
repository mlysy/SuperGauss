require(subdiff)

# log density function
logdens <- function(Z, acf1) {
  Tz <- Toeplitz(acf = acf1)
  N <- length(Z)
  IP <- colSums(Z * solve(Tz, Z))
  ldV <- determinant(Tz)
  ld <- -.5 * (IP + ldV + N * log(pi))
  ld
}

# gradient function
grad <- function(Z, acf1, dZ, dacf) {
  N <- nrow(dZ)
  p <- ncol(dZ)
  ans <- rep(NA, p)
  Tz <- Toeplitz(acf = acf1)
  T2 <- Toeplitz(n = N)
  temp <- solve(Tz, Z)
  for(ii in 1:p) {
    T2$setAcf(dacf[,ii])
    lgrad <- 2 * dZ[,ii] %*% temp - t(temp) %*% (T2 %*% temp)
    lgrad <- lgrad + Tz$traceT2(dacf[,ii])
    ans[ii] <- - lgrad / 2    
  }
  ans
}

# hessian function
hess <- function(Z, acf1, dZ, dacf, d2Z, d2acf) {
  N <- nrow(dZ)
  p <- ncol(dZ)
  ans <- matrix(NA, p, p)
  Tz <- Toeplitz(acf = acf1)
  T2 <- Toeplitz(n = N)
  temp <- solve(Tz, Z)
  for(ii in 1:p) {
    for(jj in 1:ii) {
      T2$setAcf(dacf[,jj])
      temp2 <- T2 %*% temp      
      
      T2$setAcf(dacf[,ii])
      temp1 <- T2 %*% temp
      
      lhess <- t(d2Z[,ii,jj]) %*% temp
      lhess <- lhess - t(dZ[,ii]) %*% solve(Tz, temp2)
      lhess <- lhess + t(temp1) %*% solve(Tz, temp2)
      lhess <- lhess - t(dZ[,jj]) %*% solve(Tz, temp1)
      lhess <- lhess + t(dZ[,ii]) %*% solve(Tz, dZ[,jj])
      lhess <- lhess * 2
      
      T2$setAcf(d2acf[,ii,jj])
      lhess <- lhess - t(temp) %*% (T2 %*% temp)
      
      lhess <- lhess + Tz$traceT2(d2acf[,ii,jj]) - Tz$traceT4(dacf[,ii], dacf[,jj])
      ans[ii, jj] <- -lhess / 2
    }
  }
  ans[upper.tri(ans)] <- ans[lower.tri(ans)]
  ans
}

# full gradient function w.r.t z
grad_z <- function(Z, acf1) {
  Tz <- Toeplitz(acf = acf1)
  -solve(Tz, Z)
}

upper.toep <- function(a) {
  a <- as.vector(a)
  ans <- toeplitz(a)
  ans[lower.tri(ans)] <- 0
  ans
}

# full gradient function w.r.t acf
grad_acf <- function(Z, acf1) {
  N <- length(Z)
  Tz <- Toeplitz(acf = acf1)
  Vz <- solve(Tz, Z)
  tau <- solve(Tz, c(1,rep(0, N-1)))
  ip <- upper.toep(Vz) %*% Vz
  tau2 <- c(0, tau[N:2])
  tr <- upper.toep(tau) %*% (N:1 * tau)
  tr <- tr - upper.toep(tau2) %*% (N:1 * tau2)
  tr <- tr / tau[1]
  ans <- ip - tr
  ans[1] <- ans[1] / 2
  ans
}


# test --------------------------------------------------------------------

# data generation
alpha <- .8
dT <- 1
N <- 10
p <- 3

Z <- 1:N * 1 / 10
acf1 <- fbm_acf(alpha, dT, N)
dacf <- matrix(NA, N, p)
dacf[,1] <- fbm_acf(.1, dT, N)
dacf[,2] <- fbm_acf(.2, dT, N)
dacf[,3] <- fbm_acf(.3, dT, N)
dZ <- matrix(NA, N, p)
dZ[,1] <- 1:N * .1 / 10
dZ[,2] <- 1:N * .2 / 10
dZ[,3] <- 1:N * .3 / 10
d2Z <- d2acf <- array(NA, dim = c(N, p, p))
for(ii in 1:3) {
  for(jj in 1:3) {
    d2acf[,ii,jj] <- fbm_acf((ii+jj)/10, dT, N)
    d2Z[,ii,jj] <- 1:N * (ii+jj)/100
  }
}

# log density
logdens(Z, acf1)

# gradient
grad(Z, acf1, dZ, dacf)

# hessian
c(hess(Z, acf1, dZ, dacf, d2Z, d2acf))

# gradient w.r.t. z
c(grad_z(Z, acf1))

# gradient w.r.t. acf
c(grad_acf(Z, acf1))

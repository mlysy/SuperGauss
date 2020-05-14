#--- normal with circulant variance matrix -------------------------------------

source("circulant-functions.R")

# check even fft

N <- 5
n <- floor(N/2) + 1
acf <- t2acf(N, 1:n)
acf
fft(acf)

# check acf

N <- 6
n <- floor(N/2) + 1
Ct <- Circulant$new(N = N, tacf = 1:n)
Ct$get_acf()

#--- power product -------------------------------------------------------------

# calculate y = Ct^alpha * y,

N <- 6
n <- floor(N/2) + 1
Ct <- Circulant$new(N = N, tacf = n:1)
Ct$get_acf()

N <- sample(5:20,1)
p <- sample(1:5, 1)
n <- floor(N/2) + 1
tacf <- rnorm(n)
Ct <- Circulant$new(N = N, tacf = tacf)
x <- matrix(rnorm(N*p), N, p)
y <- toeplitz(Ct$get_acf()) %*% x
y2 <- Ct$pow_prod(x)
range(y-y2)
y <- solve(toeplitz(Ct$get_acf()), x)
y2 <- Ct$pow_prod(x, pow = -1)
range(y - y2)

#--- circulant shifts and convolution ------------------------------------------

require(mvtnorm)
require(numDeriv)

N <- 5
x <- rnorm(N)

(1:N-1 + 0) %% N + 1


shift(5, 2)

N <- sample(5:20,1)
x <- rnorm(N)
y <- rnorm(N)
c1 <- sapply(1:N, function(i) crossprod(x, shift(N, i) %*% y))
c2 <- rev(convolve(x, rev(y), conj = FALSE))
c3 <- shift_mult(x, y)
range(c1 - c2)
range(c1 - c3)

#--- NormalCirculant distribution ----------------------------------------------

require(mvtnorm)
require(numDeriv)

# density
circ_ldens <- function(z, nu) {
  dmvnorm(z, sigma = toeplitz(t2acf(length(z), nu)), log = TRUE)
  ## Ct <- Circulant$new(N = length(z), tacf = nu)
  ## ld <- sum(z*Ct$solve(z)) + Ct$log_det()
  ## -.5 * ld
}

# density
N <- sample(5:20,1)
n <- floor(N/2)+1
psd <- rexp(n)
gam <- Re(FFT(t2acf(N, psd)))
nu <- gam[1:n]
z <- rnorm(N)
NCt <- NormalCirculant$new(N)
NCt$logdens(z, nu) - circ_ldens(z, nu)

# gradient calculations

# inner product term
circ_ip <- function(z, nu) {
  Ct <- Circulant$new(N = length(z), tacf = nu)
  sum(z*Ct$solve(z))
}

# trace term
circ_tr <- function(N, nu) {
  Ct <- Circulant$new(N = N, tacf = nu)
  Ct$log_det()
}

# grad wrt z
N <- sample(2:20,1)
n <- floor(N/2)+1
## nu <- Re(FFT(t2acf(N, rexp(n)*N), inverse = TRUE)[1:n]/N)
psd <- rexp(n)*N
nu <- ifft(t2acf(N, psd))[1:n]
z <- rnorm(N)
g1 <- grad(circ_ldens, x = z, nu = nu)
NCt <- NormalCirculant$new(N = N)
g2 <- NCt$grad_full(z = z, tacf = nu)$dz

y <- Ct$solve(z)
g2 <- -shift_mult(y, y)
if(N > 2) {
  eN <- (2*n) == (N+2)
  id <- 2:(n-eN)
  g2[id] <- 2*g2[id]
}
range(g1-g2[1:n])


N <- sample(2:20,1)
n <- floor(N/2)+1
## nu <- Re(FFT(t2acf(N, rexp(n)*N), inverse = TRUE)[1:n]/N)
psd <- rexp(n)*N
nu <- ifft(t2acf(N, psd))[1:n]
z <- rnorm(N)
g1 <- grad(circ_ip, x = nu, z = z)
Ct <- Circulant$new(N = N, tacf = nu)
## y <- solve(toeplitz(t2acf(N, nu)), z)
y <- Ct$solve(z)
g2 <- -shift_mult(y, y)
if(N > 2) {
  eN <- (2*n) == (N+2)
  id <- 2:(n-eN)
  g2[id] <- 2*g2[id]
}
range(g1-g2[1:n])


N <- sample(5:15, 1)
n <- floor(N/2)+1
psd <- rexp(n)*N
## gam <- Re(FFT(t2acf(N, psd), inverse = TRUE))/N
gam <- ifft(t2acf(N, psd))
nu <- gam[1:n]
V1 <- solve(toeplitz(gam))
V2 <- toeplitz(ifft(t2acf(N, 1/psd)))
range(V1 - V2)

N <- sample(5:20,1)
n <- floor(N/2)+1
psd <- rexp(n)*N
gam <- ifft(t2acf(N, psd))
nu <- gam[1:n]
g1 <- grad(circ_tr, x = nu, N = N)
## g2 <- Re(FFT(1/Re(FFT(t2acf(N, nu))), inverse = TRUE))
g2 <- ifft(1/fft(t2acf(N, nu))) * N
if(N > 2) {
  eN <- (2*n) == (N+2)
  id <- 2:(n-eN)
  g2[id] <- 2*g2[id]
}
range(g1-g2[1:n])

# put it together
N <- sample(5:20,1)
n <- floor(N/2)+1
psd <- rexp(n)*N
gam <- ifft(t2acf(N, psd))
nu <- gam[1:n]
z <- rnorm(N)
g1 <- grad(circ_ldens, x = nu, z = z)
## Ct <- Circulant$new(N = N, tacf = nu)
## y <- Ct$solve(z)
## dip <- shift_mult(y, y)
## dtr <- ifft(1/fft(t2acf(N, nu))) * N
## if(N > 2) {
##   eN <- (2*n) == (N+2)
##   id <- 2:(n-eN)
##   dip[id] <- 2*dip[id]
##   dtr[id] <- 2*dtr[id]
## }
## g2 <- .5 * (dip - dtr)[1:n]
NCt <- NormalCirculant$new(N = N)
g2 <- NCt$grad_full(z = z, tacf = nu)$dtacf
range(g1 - g2)

#--- test c++ code -------------------------------------------------------------

source("circulant-functions.R")
require(SuperGauss)

# NormalCirculant

N <- sample(5:20, 1)
Nu <- floor(N/2)+1
upsd <- rexp(Nu) * N
acf <- ifft(t2acf(N, upsd))
uacf <- acf[1:Nu]

NCt <- SuperGauss::NormalCirculant$new(N = N)
NCt2 <- NormalCirculant$new(N = N)

# log-density
z <- rnorm(N)
NCt$logdens(z = z, uacf = uacf)
NCt2$logdens(z = z, tacf = uacf)

# full gradient
grad <- NCt$grad_full(z = z, uacf = uacf)
grad2 <- NCt2$grad_full(z = z, tacf = uacf)

# wrt z
range(grad$dldz - grad2$dz)
# wrt uacf
range(grad$dldu - grad2$dtacf)

## wrt z
## dz <- -solve(toeplitz(acf), z)
## range(grad$dldz - dz)

## # check rho = ifft(fft(-dz) * fft(-rev(dz)))
## shift_mult(-dz, -dz)
## grad$dldu

## # check kappa = N * ifft(1/fft(acf))
## Re(N * ifft(1/fft(acf)))
## grad$dldu

#--- fft/ifft ------------------------------------------------------------------

N <- 10
x <- rnorm(N)
y1 <- fft(x)
y2 <- SuperGauss:::real_fft(x, inverse = FALSE)
x2 <- SuperGauss:::real_fft(y1, inverse = TRUE)

dct <- function(x) {
  n <- length(x)
  sapply(1:n-1, function(k) {
    cs <- cos(pi * (1:n-1 + 1/2)*k/n)
    2*sum(x * cs)
  })
}

lapply(setNames(2:20, 2:20), function(N) {
  n <- N %/% 2 + 1
  x <- t2acf(N, 1:n)
  ## x <- t2acf(N, rnorm(n))
  cbind(fftw = SuperGauss:::even_fft(x = c(x,x[n]))[1:n], true = fft(x)[1:n])
})


N <- sample(1:20, 1)
Nu <- N %/% 2 + 1
x <- t2acf(N, rnorm(Nu))
range(Re(fft(x))[1:Nu] - SuperGauss:::even_fft(x))
range(Re(ifft(x))[1:Nu] - SuperGauss:::even_fft(x, inverse = TRUE))

SuperGauss:::even_fft(x)
fft(x)
xp <- rep(0, 4*n)
xp <- c(1, 0, 2, 0, 3, 0, 3, 0, 2, 0)
## xp[seq(2, length(xp), by = 2)] <- x
fft(xp)

if(N %% 2 == 0) {
  cbind(fftw = SuperGauss:::even_fft(x), true = fft(x)[1:n])
} else {
  ## xp <- t2acf(N+1, c(1:n,n))
  xp <- rep(x, each = 2)
  xp[seq(2, length(xp), by = 2)] <- 0
  cbind(fftw = SuperGauss:::even_fft(xp)[1:n], true = fft(x)[1:n])
}

N <- 5
n <- N %/% 2 + 1
x <- t2acf(N, 1:n)
fft(x)
x2 <- rep(x, each = 2)
x2[seq(2, length(x2), by = 2)] <- 0
fft(x2)

EFFT_type <- setNames(0:7,
                      c("E00", "E01", "E10", "E11", "O00", "O01", "O10", "O11"))


N <- 7
n <- N %/% 2 + 1
x <- 1:n
sapply(EFFT_type, SuperGauss:::even_fft, x = x)
fft(t2acf(N, x))



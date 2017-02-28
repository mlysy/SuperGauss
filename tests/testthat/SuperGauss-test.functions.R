acf2incr.SGtest <- function(gam) {
  N <- length(gam)-1
  if(N == 1) {
    igam <- 2*(gam[1]-gam[2])
  } else {
    igam <- 2*gam[1:N] - gam[1:N+1] - gam[c(2, 1:(N-1))]
  }
  igam
}

fbm.acf.SGtest <- function(H, dT, N, incr = TRUE){
   gam <- (dT*(0:N))^(2*H)
   if(incr){
     # increments
     ans <- -1/2 * acf2incr.SGtest(gam)
   } else{
     # observations
     ans <- gam[1:N]
   }
   ans
}

d.fbm.acf.SGtest <- function(H, dT, N, incr = TRUE)

exp2.acf.SGtest <- function(lambda, dT, N, incr = TRUE) {
  # process autocorrelation
  gam <- exp(-(0:N*dT/lambda)^2)
  if(incr) {
    # increments
    ans <- acf2incr.SGtest(gam)
  } else {
    # observations
    ans <- gam[1:N]
  }
  ans
}

exp.acf.SGtest <- function(lambda, dT, N, incr = TRUE) {
  # process autocorrelation
  gam <- exp(-(0:N*dT/lambda))
  if(incr) {
    # increments
    ans <- acf2incr.SGtest(gam)
  } else {
    # observations
    ans <- gam[1:N]
  }
  ans
}

matern.acf.SGtest <- function(lambda, nu, dT, N, incr = TRUE) {
  # process autocorrelation
  tt <- sqrt(2*nu) * (0:N)*dT/lambda
  gam <- nu * log(.5 * tt) - lgamma(nu)
  gam <- 2 * exp(gam) * besselK(tt, nu)
  gam[tt == 0] <- 1
  if(incr) {
    # increments
    ans <- acf2incr.SGtest(gam)
  } else {
    # observations
    ans <- gam[1:N]
  }
  ans
}

acf.get.SGtest <- function(N, type, dT, incr = TRUE){
  lambda <- abs(rnorm(n = 1, mean = 3, sd = 1))
  H <- runif(n = 1)
  nu <- ceiling(runif(n = 1) * 4)
  if(type == "exp2"){
    exp2.acf(lambda, dT, N, incr)
  }
  if(type == "exp"){
    exp.acf(lambda, dT, N, incr)
  }
  if(type == "fbm"){
    fbm.acf(H, dT, N, incr)
  }
  if(type == "martern"){
    matern.acf(lambda, nu, dT, N, incr)
  }
  if(type == "zero"){
    rep(0, N)
  }
  if(type == "rnd"){
    rnorm(N)
  }
}

tr <- function(mat){
  if(length(mat) == 1){
    as.numeric(mat)
  }else{
    sum(diag(mat))
  }
}
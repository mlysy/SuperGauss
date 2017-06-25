# trace of matrix
.trace <- function(Mat){
  if(length(Mat) == 1){
    tr <- Mat
  }else{
    tr <- sum(diag(Mat))
  }
  tr
}

# mu should be vector of size N
.format.mu <- function(mu, N){
  if(missing(mu)){
    mu <- rep(0, N)
  }else{
    if(is.matrix(mu)){
      mu <- as.vector(mu)
    }
    if(is.vector(mu)){
      if(length(mu) != N){
        if(length(mu) == 1){
          mu <- rep(mu, N)
        }else{
          stop("mu has incompatible dimension with X")
        }
      }
    }
  }
  mu
}

# acf should be Toeplitz object of size N
.format.acf <- function(acf, N){
  if(class(acf) == "Toeplitz_Matrix"){
    # is Toeplitz
    if(ncol(acf) != N){
      stop("acf has incompatible dimension with X")
    }
  }else{
    if(is.vector(acf)){
      if(length(acf) != N){
        stop("acf has incompatible dimension with X")
      }else{
        acf <- Toeplitz(acf)
      }
    }else{
      stop("acf should be either vector or Toeplitz class")
    }
  }
  acf
}

# dmu should be matrix of size N x p, can be missing
.format.dmu <- function(dmu, N, p){
  if(missing(dmu)){
    dmu <- matrix(0, N, p)
  }else{
    if(is.vector(dmu) && p == 1){
      if(length(dmu == N)){
        dmu <- as.matrix(dmu)
      }else{
        stop("dmu has incompatible dimension with X")
      }
    }else{
      if(is.matrix(dmu)){
        if(nrow(dmu) != N || ncol(dmu) != p){
          stop("dmu has incompatible dimension with dacf")
        }
      }else{
        stop("dmu should be a matrix")
      }
    }
  }
  dmu
}

# dacf should be matrix of size N x p, cannot be missing, p decided by size of dacf
.format.dacf <- function(dacf, N){
  if(is.vector(dacf)){
    if(length(dacf == N)){
      dacf <- matrix(dacf, N, 1)
      p <- 1
    }else{
      stop("dacf has incompatible dimension with X")
    }
  }else{
    if(is.matrix(dacf)){
      if(nrow(dacf) != N){
        stop("dacf has incompatible dimension with X")
      }else{
        p <- ncol(dacf)
      }
    }else{
      stop("dacf should be a matrix")
    }
  }
  list(dacf = dacf, p = p)
}

# d2mu should be array of size N x p x p
.format.d2mu <- function(d2mu, N, p){
  if(missing(d2mu)){
    d2mu <- array(0, c(N, p, p))
  }else{
    if(is.vector(d2mu)){
      if(length(d2mu == N) && p == 1){
        d2mu <- array(d2mu, c(N, p, p))
      }else{
        stop("d2mu has incompatible dimension with X")
      }
    }else{
      if(is.array(d2mu)){
        if(!prod(dim(d2mu) == c(N, p, p))){
          stop("d2mu has incompatible dimension with dacf")
        }
      }else{
        stop("d2mu should be a array")
      }
    }
  }
  d2mu
}

# d2acf should be array of size N x p x p
.format.d2acf <- function(d2acf, N, p){
  if(is.vector(d2acf)){
    if(length(d2acf == N) && p == 1){
      d2acf <- array(d2acf, c(N, p, p))
    }else{
      stop("d2acf has incompatible dimension")
    }
  }else{
    if(is.array(d2acf)){
      if(!prod(dim(d2acf) == c(N, p, p))){
        stop("d2acf has incompatible dimension")
      }
    }else{
      stop("d2acf should be a array")
    }
  }
  d2acf
}
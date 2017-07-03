# newton-raphson for exponential covariance

# acf function
exp.acf <- function(t, lambda, sigma) sigma^2 * exp(-abs(t/lambda))
# gradient, returned as a 2-column matrix
exp.acf.grad <- function(t, lambda, sigma) {
  ea <- exp.acf(t, lambda, 1)
  cbind(abs(t)*(sigma/lambda)^2 * ea, # d_acf/d_lambda
        2*sigma * ea) # d_acf/d_sigma
}
# Hessian, returned as an array of size length(t) x 2 x 2
exp.acf.hess <- function(t, lambda, sigma) {
  ea <- exp.acf(t, lambda, 1)
  sl2 <- sigma/lambda^2
  hess <- array(NA, dim = c(length(t), 2, 2))
  hess[,1,1] <- sl2^2*(t^2 - 2*abs(t)*lambda) * ea # d2_acf/d_lambda^2
  hess[,1,2] <- 2*sl2 * abs(t) * ea # d2_acf/(d_lambda d_sigma)
  hess[,2,1] <- hess[,1,2] # d2_acf/(d_sigma d_lambda)
  hess[,2,2] <- 2 * ea # d2_acf/d_sigma^2
  hess
}

# some tests
n <- 5
tseq <- rexp(n)
lambda <- rexp(1)
sigma <- rexp(1)
a1 <- exp.acf(tseq, lambda, sigma)
ag1 <- exp.acf.grad(tseq, lambda, sigma)
ah1 <- exp.acf.hess(tseq, lambda, sigma)
ii <- sample(n, 1)
ag2 <- grad(func = function(ls) exp.acf(tseq[ii], ls[1], ls[2]),
            x = c(lambda, sigma))
ah2 <- hessian(func = function(ls) exp.acf(tseq[ii], ls[1], ls[2]),
               x = c(lambda, sigma))
ag1[ii,]-ag2
ah1[ii,,]-ah2

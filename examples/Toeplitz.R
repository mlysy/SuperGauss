# construct a Toeplitz matrix
acf <- exp(-(1:5))
Tz <- Toeplitz$new(acf = acf)
# alternatively, can allocate space first
Tz <- Toeplitz$new(N = length(acf))
Tz$set_acf(acf = acf)

# basic methods
Tz$get_acf() # extract the acf
dim(Tz) # == c(nrow(Tz), ncol(Tz))
Tz # print method

# linear algebra methods
X <- matrix(rnorm(10), 5, 2)
Tz %*% X
t(X) %*% Tz
solve(Tz, X)
determinant(Tz) # log-determinant

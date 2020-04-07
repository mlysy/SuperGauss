N <- 10
acf <- exp(-(1:N - 1)/N)
rnormtz(n = 3, acf = acf)

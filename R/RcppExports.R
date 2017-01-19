# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

toeplitzXZ <- function(X, acf) {
    .Call('SuperGauss_DurbinLevinson_XZ', PACKAGE = 'SuperGauss', X, acf)
}

toeplitzZX <- function(Z, acf) {
    .Call('SuperGauss_DurbinLevinson_ZX', PACKAGE = 'SuperGauss', Z, acf)
}

DurbinLevinsonEigen <- function(X, Y, acf, calcMode = 1L) {
    .Call('SuperGauss_DurbinLevinson_Eigen', PACKAGE = 'SuperGauss', X, Y, acf, calcMode)
}

DurbinLevinsonBase <- function(X, Y, acf, calcMode = 1L) {
    .Call('SuperGauss_DurbinLevinson_Base', PACKAGE = 'SuperGauss', X, Y, acf, calcMode)
}


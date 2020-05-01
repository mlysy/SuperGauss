#--- R6 tests ------------------------------------------------------------------

require(SuperGauss)
require(R6)
require(Rcpp)

sourceCpp("r6_test.cpp")

foo <- R6Class(
  classname = "foo",
  public = list(
    vf_ = NULL,
    g2k_ = NULL,
    gn_ = NULL,
    initialize = function(N) {
      ## self$vf_ <- VectorFFT_ctor(N)
      ## self$g2k_ <- GSchur2K_ctor(N)
      self$gn_ <- GSchurN_ctor(N)
    }
  )
)

## bar <- foo$new(4e5)
N <- 4e5
bar <- GSchurN_ctor(N)
rm(bar); gc()

replicate(50, {
  ## bar <- foo$new(4e5)
  bar <- GSchurN_ctor(N)
  rm(bar); gc()
})

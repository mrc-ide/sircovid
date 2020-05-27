## Compile the shared library with
##
## R CMD SHLIB unuran_binom.c
dyn.load("unuran_binom.so")

test_unran_binom_r <- function(n, p) {
  .Call("test_unuran_binom_r", as.integer(n), as.numeric(p), PACKAGE = "unuran_binom")
}

dyn.load("binom.so")

test_binom_r <- function(n, p) {
  .Call("test_binom_r", as.integer(n), as.numeric(p), PACKAGE = "binom")
}

m <- 1e5
n <- rpois(m, 50)
p <- runif(m)

bench::mark(
  test_unran_binom_r(n, p),
  test_binom_r(n, p),
  rbinom(m, n, p),
  check = FALSE)

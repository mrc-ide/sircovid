## Compile the shared library with
##
## R CMD SHLIB binom.c
dyn.load("binom.so")

test_binom_r <- function(n, p) {
  .Call("test_binom_r", as.integer(n), as.numeric(p), PACKAGE = "binom")
}

m <- 1e5
n <- rpois(m, 50)
p <- runif(m)
test_binom_r(n, p)

bench::mark(
  test_binom_r(n, p),
  rbinom(m, n, p),
  check = FALSE)

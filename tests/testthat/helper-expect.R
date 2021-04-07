expect_rounded_lte <- function(x, y, digits = 10, tol = 1e-5) {
 expect_true(all(round(x, digits) - round(y, digits) <= tol))
}

expect_rounded_lt <- function(x, y, digits = 10, tol = 1e-5) {
 expect_true(all(round(x, digits) - round(y, digits) < tol))
}

expect_rounded_gte <- function(x, y, digits = 10, tol = -1e-5) {
 expect_true(all(round(x, digits) - round(y, digits) >= tol))
}

expect_rounded_gt <- function(x, y, digits = 10, tol = -1e-5) {
 expect_true(all(round(x, digits) - round(y, digits) > tol))
}

expect_vector_equal <- function(x, y, digits = 100, tol = 0) {
 if (is.numeric(x) && is.numeric(y)) {
   expect_true(all(abs(round(x, digits) - round(y, digits)) <= tol),
              sprintf("\nNot all %s equal to %s (tol %s)",
              deparse(substitute(x)), deparse(substitute(y)), tol))
 } else {
   expect_true(all(x == y))
 }
}


expect_vector_lte <- function(x, y, digits = 100, tol = 0) {
 expect_true(all(round(x, digits) - round(y, digits) <= tol),
              sprintf("\nNot all %s less than or equal to %s (tol %s)",
              deparse(substitute(x)), deparse(substitute(y)), tol))
}


expect_vector_lt <- function(x, y, digits = 100, tol = 0) {
 expect_true(all(round(x, digits) - round(y, digits) < tol),
              sprintf("\nNot all %s less than %s (tol %s)",
              deparse(substitute(x)), deparse(substitute(y)), tol))
}


expect_vector_gte <- function(x, y, digits = 100, tol = -0) {
 expect_true(all(round(x, digits) - round(y, digits) >= tol),
              sprintf("\nNot all %s greater than or equal to %s (tol %s)",
              deparse(substitute(x)), deparse(substitute(y)), tol))
}


expect_vector_gt <- function(x, y, digits = 100, tol = -0) {
 expect_true(all(round(x, digits) - round(y, digits) > tol),
              sprintf("\nNot all %s greater than %s (tol %s)",
              deparse(substitute(x)), deparse(substitute(y)), tol))
}

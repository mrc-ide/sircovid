#------------------------------------------------
# x is what
#' @noRd
assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name,
                 paste(what, collapse = " / ")), call. = FALSE)
  }
}

#------------------------------------------------
# x is positive integer (with or without zero allowed)
#' @noRd
assert_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_int(x, name = name)
  assert_pos(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}
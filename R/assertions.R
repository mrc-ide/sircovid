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
# x is numeric
#' @noRd
assert_numeric <- function(x, message = "%s must be numeric", name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is integer
#' @noRd
assert_int <- function(x, message = "%s must be integer valued", name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (!isTRUE(all.equal(x, as.integer(x), check.attributes = FALSE))) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}


#------------------------------------------------
# x is positive (with or without zero allowed)
#' @noRd
assert_pos <- function(x, zero_allowed = TRUE, message1 = "%s must be greater than or equal to zero", message2 = "%s must be greater than zero", name = deparse(substitute(x))) {
  assert_numeric(x, name = name)
  if (zero_allowed) {
    if (!all(x>=0)) {
      stop(sprintf(message1, name), call. = FALSE)
    }
  } else {
    if (!all(x>0)) {
      stop(sprintf(message2, name), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# x is positive integer (with or without zero allowed)
#' @noRd
assert_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_int(x, name = name)
  assert_pos(x, zero_allowed = zero_allowed, name = name)
  return(TRUE)
}
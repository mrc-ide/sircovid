assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name, paste(what, collapse = " / ")),
         call. = FALSE)
  }
  invisible(x)
}


assert_increasing <- function(x, strict = TRUE, name = deparse(substitute(x))) {
  dx <- diff(x)
  if (strict) {
    if (any(dx <= 0)) {
      stop(sprintf("'%s' must be strictly increasing", name))
    }
  } else {
    if (any(dx < 0)) {
      stop(sprintf("'%s' must be increasing (or equal)", name))
    }
  }
  invisible(x)
}

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


assert_integer <- function(x, name = deparse(substitute(x))) {
  if (!all(is_integer(x))) {
    stop(sprintf("'%s' must be an integer", name), call. = FALSE)
  }
  invisible(x)
}


assert_date_string <- function(x, name = deparse(substitute(x))) {
  if (inherits(x, "Date")) {
    x <- as.character(x)
  }
  if (!all(grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", x))) {
    stop(sprintf("Expected ISO dates or R dates for '%s' - please convert",
                 name))
  }
  invisible(x)
}


assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}


assert_non_negative <- function(x, name = deparse(substitute(x))) {
  if (any(x < 0)) {
    stop(sprintf("'%s' must have only non-negative values", name),
         call. = FALSE)
  }
  invisible(x)
}


assert_positive <- function(x, name = deparse(substitute(x))) {
  if (any(x <= 0)) {
    stop(sprintf("'%s' must have only positive values", name),
         call. = FALSE)
  }
  invisible(x)
}


assert_proportion <- function(x, name = deparse(substitute(x))) {
  if (any(x < 0 | x > 1)) {
    stop(sprintf("'%s' must lie in [0, 1]", name),
         call. = FALSE)
  }
  invisible(x)
}


assert_relatives <- function(x, name = deparse(substitute(x))) {
  if (x[[1]] != 1) {
    stop(sprintf("'%s[[1]]' must be 1", name), call. = FALSE)
  }
  assert_non_negative(x, name)
  invisible(x)
}


assert_logical <- function(x, name = deparse(substitute(x))) {
  if (!(all(x %in% 0:1) || is.logical(x))) {
    stop(sprintf("All '%s' must be logical or in {0, 1}", name), call. = FALSE)
  }

  invisible(as.integer(x))
}

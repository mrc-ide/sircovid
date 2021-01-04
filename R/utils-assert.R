assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name, paste(what, collapse = " / ")),
         call. = FALSE)
  }
  invisible(x)
}


assert_increasing <- function(x, name = deparse(substitute(x))) {
  if (any(diff(x) <= 0)) {
    stop(sprintf("'%s' must be strictly increasing", name))
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

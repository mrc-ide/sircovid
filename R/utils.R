`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}


squote <- function(x) {
  sprintf("'%s'", x)
}


verify_names <- function(x, required = NULL, optional = NULL,
                         allow_extra = FALSE,
                         name = deparse(substitute(x))) {
  nms <- names(x)
  if (anyDuplicated(nms)) {
    dups <- unique(nms[duplicated(nms)])
    stop(sprintf("Duplicate element names in '%s': %s",
                 dups, paste(squote(x), collapse = ", ")))
  }
  if (!allow_extra) {
    extra <- setdiff(nms, c(required, optional))
    if (length(extra) > 0) {
      stop(sprintf("Extra elements in '%s': %s",
                   name, paste(squote(x), collapse = ", ")))
    }
  }
  msg <- setdiff(required, nms)
  if (length(msg) > 0) {
    stop(sprintf("Elements missing from '%s': %s",
                 name, paste(squote(x), collapse = ", ")))
  }
  invisible(x)
}


vlapply <- function(fun, x, ...) {
  vapply(fun, x, TRUE, ...)
}


vnapply <- function(fun, x, ...) {
  vapply(fun, x, 1.0, ...)
}


sircovid_file <- function(...) {
  system.file(..., package = "sircovid2", mustWork = TRUE)
}


read_csv <- function(...) {
  utils::read.csv(..., check.names = FALSE, stringsAsFactors = FALSE)
}


data_frame <- function(...) {
  data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
}


is_integer <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}


rename <- function(x, from, to, name = deparse(substitute(x))) {
  verify_names(x, required = from, allow_extra = TRUE, name = name)
  i <- match(from, names(x))
  names(x)[i] <- to
  x
}


is_date <- function(x) {
  inherits(date, "Date")
}


as_date <- function(date) {
  if (is_date(date)) {
    return(date)
  }
  if (!all(grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", date))) {
    stop("Expected ISO dates or R dates - please convert")
  }
  as.Date(date)
}

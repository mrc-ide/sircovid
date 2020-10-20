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
                 name, paste(squote(dups), collapse = ", ")))
  }
  if (!allow_extra) {
    extra <- setdiff(nms, c(required, optional))
    if (length(extra) > 0) {
      stop(sprintf("Extra elements in '%s': %s",
                   name, paste(squote(extra), collapse = ", ")))
    }
  }
  msg <- setdiff(required, nms)
  if (length(msg) > 0) {
    stop(sprintf("Elements missing from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  invisible(x)
}


sircovid_file <- function(...) {
  system.file(..., package = "sircovid", mustWork = TRUE)
}


read_csv <- function(...) {
  utils::read.csv(..., check.names = FALSE, stringsAsFactors = FALSE)
}


write_csv <- function(...) {
  utils::write.csv(..., row.names = FALSE)
}


data_frame <- function(...) {
  data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
}


is_integer <- function(x, tol = sqrt(.Machine$double.eps)) {
  abs(x - round(x)) < tol
}


rename <- function(x, from, to, name = deparse(substitute(x))) {
  verify_names(x, required = from, allow_extra = TRUE, name = name)
  i <- match(from, names(x))
  names(x)[i] <- to
  x
}


set_names <- function(x, nms) {
  names(x) <- nms
  x
}


vnapply <- function(x, fun, ...) {
  vapply(x, fun, numeric(1), ...)
}

check_rel_susceptibility <- function(rel_susceptibility) {
  between_zero_and_one <- all(rel_susceptibility >= 0 & rel_susceptibility <= 1)
  first_group_is_reference <- rel_susceptibility[1] == 1
  if(length(rel_susceptibility) > 1) {
    # relative susceptibility is lowest in vaccinated group and then increases
    waning_immunity <- all(diff(rel_susceptibility[-1]) > 0)
  } else waning_immunity <- TRUE
  message <- paste("rel_susceptibility does not have the expected format i.e.:",
                   "- all values between 0 and 1, ",
                   "- first value is 1,",
                   "- and values thereafter are increasing.",
                   sep = "\n")
  if(!(between_zero_and_one & first_group_is_reference & waning_immunity))
    stop(message)
}
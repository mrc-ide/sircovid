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

build_rel_susceptibility <- function(rel_susceptibility, N_age) {
  if (is.matrix(rel_susceptibility)) {
    if (nrow(rel_susceptibility) != N_age) {
      stop("'rel_susceptibility' should have as many rows as age groups")
    }
    for (i in seq_len(nrow(rel_susceptibility))) {
      check_rel_susceptibility(rel_susceptibility[i, ])
    }
    mat_rel_susceptibility <- rel_susceptibility
  } else { # create matrix by repeating rel_susceptibility for each age group
    mat_rel_susceptibility <-
      matrix(rep(rel_susceptibility, each = N_age), nrow = N_age)
  }
  mat_rel_susceptibility
}

check_rel_susceptibility <- function(rel_susceptibility) {
  if (length(rel_susceptibility) == 0) {
    stop("At least one value required for 'rel_susceptibility'")
  }
  if (any(rel_susceptibility < 0 | rel_susceptibility > 1)) {
    stop("All values of 'rel_susceptibility' must lie in [0, 1]")
  }
  if (rel_susceptibility[[1]] != 1) {
    stop("First value of 'rel_susceptibility' must be 1")
  }
}

## TO DO: add checks that 

# vaccination_rate is a vector of length N_age
# (eventually needs to make it possible to vary over time)

# vaccine_progression_rate is a matrix 
# with N_age rows and N_vacc_classes - 2 columns


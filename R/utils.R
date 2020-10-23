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

build_vaccination_rate <- function(vaccination_rate, N_age, N_vacc_classes) {
  if (N_vacc_classes == 1) { # no vaccination
    vect_vaccination_rate <- rep(0, N_age)
  } else {
    if (length(vaccination_rate) > 1) {
      if (length(vaccination_rate) != N_age) {
        stop("'vaccination_rate' should have as many elements as age groups")
      }
      if (any(vaccination_rate < 0)) {
        stop("'vaccination_rate' must have only non-negative values")
      }
      vect_vaccination_rate <- vaccination_rate
    } else { # create vector by repeating vaccination_rate for each age group
      vect_vaccination_rate <- rep(vaccination_rate, each = N_age)
    }
  }
  vect_vaccination_rate
}

build_vaccine_progression_rate <- function(vaccine_progression_rate,
                                           N_age, N_vacc_classes) {
  if (N_vacc_classes <= 2) { # no vaccine progression
    mat_vaccine_progression_rate <- matrix(0, N_age, 1)
  }
  else {
    if (is.matrix(vaccine_progression_rate)) {
      if (nrow(vaccine_progression_rate) != N_age) {
        stop("'vaccine_progression_rate' must have as many rows as age groups")
      }
      if (ncol(vaccine_progression_rate) != N_vacc_classes - 2) {
        stop(
          "'vaccine_progression_rate' must have 'N_vacc_classes - 2' columns")
      }
      if (any(vaccine_progression_rate < 0)) {
        stop("'vaccine_progression_rate' must have only non-negative values")
      }
      mat_vaccine_progression_rate <- vaccine_progression_rate
    } else { # vaccine_progression_rate is a vector of length N_vacc_classes - 2
      if (!is.vector(vaccine_progression_rate) ||
        length(vaccine_progression_rate) != N_vacc_classes - 2) {
        msg1 <- "'vaccine_progression_rate' must be either:"
        msg2 <- "a vector of length 'N_vacc_classes - 2'"
        msg3 <- "or a matrix with 'N_age' rows and 'N_vacc_classes - 2' columns"
        stop(paste(msg1, msg2, msg3))
      }
      if (any(vaccine_progression_rate < 0)) {
        stop("'vaccine_progression_rate' must have only non-negative values")
      }
      # create matrix by repeating vaccine_progression_rate for each age group
      mat_vaccine_progression_rate <-
        matrix(rep(vaccine_progression_rate, each = N_age), nrow = N_age)
    }
  }
  mat_vaccine_progression_rate
}

build_rel_param <- function(rel_param, name_param) {
  n_groups <- carehomes_n_groups()
  if (is.matrix(rel_param)) {
    if (nrow(rel_param) != n_groups) {
      stop(paste(name_param, "should have as many rows as age groups"))
    }
    for (i in seq_len(nrow(rel_param))) {
      check_rel_param(rel_param[i, ], name_param)
    }
    mat_rel_param <- rel_param
  } else { # create matrix by repeating rel_param for each age group
    mat_rel_param <-
      matrix(rep(rel_param, each = n_groups), nrow = n_groups)
  }
  mat_rel_param
}


check_rel_param <- function(rel_param, name_param) {
  if (length(rel_param) == 0) {
    stop(paste("At least one value required for", name_param))
  }
  if (any(rel_param < 0 | rel_param > 1)) {
    stop(paste("All values of", name_param, "must lie in [0, 1]"))
  }
  if (rel_param[[1]] != 1) {
    stop(paste("First value of", name_param, "must be 1"))
  }
}


build_vaccine_progression_rate <- function(vaccine_progression_rate,
                                           n_vacc_classes) {
  n_groups <- carehomes_n_groups()
  # if NULL, set vaccine_progression_rate to 0
  if (is.null(vaccine_progression_rate)) {
    mat_vaccine_progression_rate <- matrix(0, n_groups, n_vacc_classes)
  } else {
    if (is.matrix(vaccine_progression_rate)) {
      if (nrow(vaccine_progression_rate) != n_groups) {
        stop(
          "'vaccine_progression_rate' must have as many rows as age groups")
      }
      if (ncol(vaccine_progression_rate) != n_vacc_classes) {
        stop(
          "'vaccine_progression_rate' must have 'n_vacc_classes' columns")
      }
      if (any(vaccine_progression_rate < 0)) {
        stop("'vaccine_progression_rate' must have only non-negative values")
      }
      mat_vaccine_progression_rate <- vaccine_progression_rate
    } else { # vaccine_progression_rate vector of length n_vacc_classes
      if (!is.vector(vaccine_progression_rate) ||
          length(vaccine_progression_rate) != n_vacc_classes) {
        m1 <- "'vaccine_progression_rate' must be either:"
        m2 <- "a vector of length 'n_vacc_classes' or"
        m3 <- "a matrix with 'n_groups' rows and 'n_vacc_classes' columns"
        stop(paste(m1, m2, m3))
      }
      if (any(vaccine_progression_rate < 0)) {
        stop("'vaccine_progression_rate' must have only non-negative values")
      }
      # create matrix by repeating vaccine_progression_rate for each age group
      mat_vaccine_progression_rate <-
        matrix(rep(vaccine_progression_rate, each = n_groups), nrow = n_groups)
    }
  }
  if (n_vacc_classes == 1 & !all(mat_vaccine_progression_rate == 0)) {
    msg1 <- "When 'n_vacc_classes' is 1,"
    msg2 <- "'vaccine_progression_rate' should only contain zeros"
    stop(paste(msg1, msg2))
  }
  mat_vaccine_progression_rate
}

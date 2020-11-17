build_rel_param <- function(rel_param, n_vacc_classes, name_param) {
  n_groups <- carehomes_n_groups()
  if (length(rel_param) == 1) {
    mat_rel_param <- matrix(rel_param, n_groups, n_vacc_classes)
  } else if (is.matrix(rel_param)) {
    if (nrow(rel_param) != n_groups) {
      stop(paste(name_param, "should have as many rows as age groups"))
    }
    mat_rel_param <- rel_param
  } else { # create matrix by repeating rel_param for each age group
    mat_rel_param <-
      matrix(rep(rel_param, each = n_groups), nrow = n_groups)
  }
  check_rel_param(mat_rel_param, name_param)
  mat_rel_param
}


check_rel_param <- function(rel_param, name_param) {
  if (length(rel_param) == 0) {
    stop(paste("At least one value required for", name_param))
  }
  if (any(rel_param < 0 | rel_param > 1)) {
    stop(paste("All values of", name_param, "must lie in [0, 1]"))
  }
  if (!all(rel_param[, 1] == 1)) {
    stop(paste("First value of", name_param, "must be 1"))
  }
}


build_time_varying_vaccine_progression_rate <- function(date, value, dt,
                                                        n_vacc_classes) {
  if (is.null(date)) {
    if (is.list(value)) {
      stop(paste("As 'date' is NULL, expected a vector or matrix for",
                 "vaccine_progression_rate"))
    }
    ret <- build_vaccine_progression_rate(value, n_vacc_classes)
    return(array(ret, c(1, dim(ret))))
  }
  if (date[[1]] != 0) {
    stop("Unlike beta, the first date has no effect and must be zero")
  }
  value <- lapply(value, build_vaccine_progression_rate, n_vacc_classes)
  dat <- sircovid_parameters_time_varying(date, value, dt,
                                          "vaccine_progression_rate")

  len <- c(dim(dat$value[[1]]), length(dat$value))
  value <- aperm(array(unlist(dat$value), len), c(3, 1, 2))

  date_vec <- seq(0, max(date), by = dt)
  value[findInterval(date_vec, date), , , drop = FALSE]
}


build_vaccine_progression_rate <- function(vaccine_progression_rate,
                                           n_vacc_classes) {
  n_groups <- carehomes_n_groups()
  # if NULL, set vaccine_progression_rate to 0
  if (is.null(vaccine_progression_rate)) {
    mat_vaccine_progression_rate <- matrix(0, n_groups, n_vacc_classes)
  } else if (is.matrix(vaccine_progression_rate)) {
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
  if (n_vacc_classes == 1 & !all(mat_vaccine_progression_rate == 0)) {
    msg1 <- "When 'n_vacc_classes' is 1,"
    msg2 <- "'vaccine_progression_rate' should only contain zeros"
    stop(paste(msg1, msg2))
  }
  mat_vaccine_progression_rate
}


get_n_candidates_vaccine_progression <- function(carehomes_dust_model_output,
                                                 mod) {
  y <- mod$transform_variables(drop(carehomes_dust_model_output))
  unvacc_stage <- 1
  S_candidates <- y$S[, unvacc_stage, ]
  ## TODO: write the two next lines in more generic a way which does not assume
  ## that dim(y$E)[2] > 1 and that dim(y$I_asympt)[2] = 1
  E_candidates <- apply(y$E[, , unvacc_stage, ], c(1, 3), sum)
  I_asympt_candidates <- drop(y$I_asympt[, , unvacc_stage, ])
  R_candidates <- drop(y$R[, unvacc_stage, ])

  S_candidates + E_candidates + I_asympt_candidates + R_candidates
}
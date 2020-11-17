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


build_time_varying_vaccine_progression_rate <- function(
  vaccine_progression_rate_date,
  vaccine_progression_rate_value, n_groups, n_vacc_classes) {
  
  ## assumes vaccine_progression_rate_value is an array with dimensions
  ## n_groups, n_vacc_classes, length(vaccine_progression_rate_date)
  ## this means we need to have called build_vaccine_progression_rate before
  
  if (is.null(vaccine_progression_rate_value))
  {
    if (!is.null(vaccine_progression_rate_date) && length(vaccine_progression_rate_date) > 1) {
      msg1 <- "if 'vaccine_progression_rate_value' is NULL"
      msg2 <- "vaccine_progression_rate_date should contain at most one date"
      stop(paste(msg1, msg2))
    }
    vaccine_progression_rate_value <- array(0, c(n_groups, n_vacc_classes, 1))
  }
  
  dt <- 0.25
  i <- 1
  j <- 1
  # to get the correct dimension for the output
  tmp <- sircovid_parameters_time_varying(
    vaccine_progression_rate_date,
    vaccine_progression_rate_value[i, j, ] %||% 0, dt)
  
  vaccine_progression_rate_step <- 
    array(NA, c(dim(vaccine_progression_rate_value)[1], 
                dim(vaccine_progression_rate_value)[2], 
                length(tmp)))
  
  for (i in seq_len(dim(vaccine_progression_rate_value)[1])) {
    for (j in seq_len(dim(vaccine_progression_rate_value)[2])) {
      value <- sapply(seq_len(dim(vaccine_progression_rate_value)[3]),
                      function(k) vaccine_progression_rate_value[i, j, k])
      vaccine_progression_rate_step[i, j, ] <- 
        sircovid_parameters_time_varying(vaccine_progression_rate_date, 
                                         value %||% 0, dt)
    }
  }
  
  vaccine_progression_rate_step
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

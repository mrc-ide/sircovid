##' Remap state that was run with one stratum to support
##' vaccination. This is useful where fitting was done pre-vaccination
##' and the model state needs to be expanded to support vaccination
##' for simulations.
##'
##' @title Remap state to include vaccination
##'
##' @param state_orig An original state matrix (rows representing
##'   state, columns representing particles or samples). This must be
##'   the complete model state.
##'
##' @param info_orig The results of `info()` from the model object
##'   without vaccination
##'
##' @param info_vacc The results of `info()` from the model object
##'   with vaccination
##'
##' @return A 3d array with more rows than `state`.
##' @export
vaccination_remap_state <- function(state_orig, info_orig, info_vacc) {
  state_vacc <- matrix(0.0, info_vacc$len, ncol(state_orig))

  extra <- setdiff(names(info_orig$index), names(info_vacc$index))
  if (length(extra) > 0) {
    stop(sprintf("Can't downgrade state (previously had variables %s)",
                 paste(squote(extra), collapse = ", ")))
  }
  msg <- setdiff(names(info_vacc$index), names(info_orig$index))
  ## We can tolerate any vaccine-related variable that can start
  ## zero'd. This will be the case through brief windows of upgrading
  ## sircovid only (e.g., between sircovid 0.7.2 and 0.8.0)
  allowed <- "cum_n_vaccinated"
  err <- setdiff(msg, allowed)
  if (length(err) > 0) {
    stop(sprintf("Can't remap state (can't add variables %s)",
                 paste(squote(err), collapse = ", ")))
  }

  for (nm in names(info_orig$index)) {
    i_orig <- info_orig$index[[nm]]
    i_vacc <- info_vacc$index[[nm]]
    if (length(i_orig) == length(i_vacc)) {
      state_vacc[i_vacc, ] <- state_orig[i_orig, ]
    } else {
      d_orig <- info_orig$dim[[nm]]
      d_vacc <- info_vacc$dim[[nm]]
      nd <- length(d_orig)
      j <- seq_len(prod(d_orig[-nd]))
      state_vacc[i_vacc[j], ] <- state_orig[i_orig, ]
    }
  }

  state_vacc
}


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
  if (!all(mat_vaccine_progression_rate[, 1] == 0)) {
    stop("The first column of 'vaccine_progression_rate' must be zero")
  }
  mat_vaccine_progression_rate
}

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
vaccine_remap_state <- function(state_orig, info_orig, info_vacc) {
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
  allowed <- c("cum_n_vaccinated", "D",
               "cum_infections_disag", "diagnoses_admitted")
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


build_rel_param <- function(rel_param, n_strains, n_vacc_classes, name_param) {
  n_groups <- carehomes_n_groups()
  if (length(rel_param) == 1) {
    mat_rel_param <- array(rel_param, c(n_groups, n_strains, n_vacc_classes))
  } else if (is.array(rel_param)) {
    if (length(dim(rel_param)) != 3) {
      stop(paste(name_param, "should be a three dimensional array with",
                 "dimensions: age groups, strains, vaccine classes"))
    }
    if (nrow(rel_param) != n_groups) {
      stop(paste(name_param, "should have as many rows as age groups"))
    }
    if (ncol(rel_param) != n_strains) {
      stop(paste(name_param, "should have as many columns as strains"))
    }
    if (dim(rel_param)[3] != n_vacc_classes) {
      stop(paste(name_param,
                 "should have number of vaccine classes as 3rd dimension"))
    }
    mat_rel_param <- rel_param
  } else { # create array by repeating rel_param for each age group and strain
    mat_rel_param <-
      array(rep(rel_param, each = n_groups * n_strains),
            dim = c(n_groups, n_strains, n_vacc_classes))
  }
  check_rel_param(mat_rel_param, name_param)
  mat_rel_param
}


check_rel_param <- function(rel_param, name_param) {
  if (length(rel_param) == 0) {
    stop(paste("At least one value required for", name_param))
  }
  if (ncol(rel_param) == 1 && !all(rel_param[, , 1] == 1)) {
    stop(paste("First value of", name_param, "must be 1"))
  } else if (ncol(rel_param) > 1 &&
             !all(rel_param[, 1, 1, drop = FALSE] == 1)) {
    stop(paste("First value of", name_param,
               "must be 1 for first infection with strain one"))
  }
}


build_vaccine_progression_rate <- function(vaccine_progression_rate,
                                           n_vacc_classes,
                                           index_dose) {
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
  for (i in seq_along(index_dose)) {
    j <- index_dose[[i]]
    if (!all(mat_vaccine_progression_rate[, j] == 0)) {
      stop(sprintf(
        "Column %d of 'vaccine_progression_rate' must be zero (dose %d)",
        j, i))
    }
  }
  mat_vaccine_progression_rate
}


##' Compute vaccination priority following JCVI ordering.
##'
##' https://tinyurl.com/8uwtatvm
##'
##' Assmuing independance between job (e.g. HCW) and clinical condition
##'
##' But assuming one is either counted as "clinically extremely
##' vulnerable" or "with underlying health conditions" but not both
##'
##' The JCVI priority groups, in decending order are:
##'
##'  1. residents in a care home for older adults and their carers
##'  2. all those 80 years of age and over and frontline health and social
##'     care workers
##'  3. all those 75 years of age and over
##'  4. all those 70 years of age and over and clinically extremely
##'     vulnerable individuals
##'  5. all those 65 years of age and over
##'  6. all individuals aged 16 years to 64 years with underlying health
##'     conditions which put them at higher risk of serious disease and
##'     mortality
##'  7. all those 60 years of age and over
##'  8. all those 55 years of age and over
##'  9. all those 50 years of age and over
##' 10. all those 40-49 years of age and over
##' 11. all those 30-39 years of age and over
##' 12. all those 18-29 years of age and over
##' @title Compute vaccination order
##'
##' @param uptake A vector of length 19 with fractional uptake per
##'   group. If a single number is given it is shared across all
##'   groups (note that this includes under-18s)
##'
##' @param prop_hcw Assumed fraction of healthcare workers in each
##'   group (length 19) - if `NULL` we use a default that is a guess
##'   with hopefully the right general shape.
##'
##' @param prop_very_vulnerable Assumed fraction "very vulnerable" in
##'   each group (length 19) - if `NULL` we use a default that is a
##'   guess with hopefully the right general shape.
##'
##' @param prop_underlying_condition Assumed fraction "underlying
##'   condition" in each group (length 19) - if `NULL` we use a
##'   default that is a guess with hopefully the right general shape.
##'
##' @return A matrix with n_groups rows (19) and columns representing
##'   priority groups, and element (i, j) is the proportion in group i
##'   who should be vaccinated as part of priority group j, accounting
##'   for uptake so that the sum over the rows corresponds to the
##'   total fractional uptake in that group. For
##'   `vaccine_priority_population`, the total number of
##'   individuals replaces the proportion (based on the demography
##'   used by sircovid).
##'
##' @rdname vaccine_priority
##' @export
vaccine_priority_proportion <- function(uptake,
                                        prop_hcw = NULL,
                                        prop_very_vulnerable = NULL,
                                        prop_underlying_condition = NULL) {
  n_groups <- carehomes_n_groups()

  if (is.null(uptake)) {
    uptake <- rep(1, n_groups)
  } else {
    uptake <- recycle(uptake, n_groups)
  }

  prop_hcw <- prop_hcw %||%
    c(rep(0, 4), rep(0.1, 10), rep(0, 5))
  prop_very_vulnerable <- prop_very_vulnerable %||%
    c(rep(0, 4), rep(0.05, 5), rep(0.1, 5), rep(0.15, 5))
  prop_underlying_condition <- prop_underlying_condition %||%
    c(rep(0, 4), rep(0.05, 5), rep(0.1, 5), rep(0.15, 5))

  n_priority_groups <- 12
  p <- matrix(0, n_groups, n_priority_groups)

  ## Aged base priority list
  jcvi_priority <- list(
    ## the age groups targeted in each priority group (see comments above)
    18:19, 17, 16, 15, 14, NULL, 13, 12, 11, 9:10, 7:8, 1:6)

  ## 1. Start with non aged based priority:
  ## helper function
  add_prop_to_vacc <- function(j, idx, prop_to_vaccinate, p) {
    p[idx, j] <- prop_to_vaccinate[idx]
    p
  }

  ## Group 2 includes frontline health and social care workers
  p <- add_prop_to_vacc(j = 2,
                        idx = seq_len(jcvi_priority[[2]] - 1),
                        prop_to_vaccinate = prop_hcw,
                        p)

  ## Group 4 includes clinically extremely vulnerable individuals
  p <- add_prop_to_vacc(j = 4,
                        idx = seq_len(jcvi_priority[[4]] - 1),
                        prop_to_vaccinate = (1 - prop_hcw) *
                          prop_very_vulnerable,
                        p)

  ## Group 6 includes all individuals aged 16 years to 64 years with
  ## underlying health conditions which put them at higher risk of
  ## serious disease and mortality
  p <- add_prop_to_vacc(j = 6,
                        idx = 4:13,
                        prop_to_vaccinate = (1 - prop_hcw) *
                          prop_underlying_condition,
                        p)

  ## 2. Add aged base priority
  for (j in seq_along(jcvi_priority)) {
    if (!is.null(jcvi_priority[[j]])) {
      ## discount those already vaccinated as part of non age based priority
      p[jcvi_priority[[j]], j] <- 1 - rowSums(p)[jcvi_priority[[j]]]
    }
  }

  ## 3. Account for uptake
  uptake_mat <- matrix(rep(uptake, n_priority_groups),
                       nrow = n_groups)

  p * uptake_mat
}


##' @param region Region to use to get total population numbers
##'
##' @rdname vaccine_priority
##' @export
vaccine_priority_population <- function(region,
                                        uptake,
                                        prop_hcw = NULL,
                                        prop_very_vulnerable = NULL,
                                        prop_underlying_condition = NULL) {
  p <- vaccine_priority_proportion(uptake,
                                   prop_hcw,
                                   prop_very_vulnerable,
                                   prop_underlying_condition)
  ## TODO: it would be nice to make this easier
  pop <- carehomes_parameters(1, region)$N_tot
  pop_mat <- matrix(rep(pop, ncol(p)), nrow = nrow(p))
  round(p * pop_mat)
}


##' Create future vaccination schedule from a projected number of
##' daily doses.
##'
##' @title Create vaccination schedule
##'
##' @param start Either a [sircovid_date] object corresponding to the
##'   first date in daily_doses_value, or a [vaccine_schedule] object
##'   corresponding to previously carried out vaccination.
##'
##' @param daily_doses_value A vector of doses per day.
##'
##' @param mean_days_between_doses Assumed mean days between doses one
##'   and two
##'
##' @param priority_population Output from
##'   [vaccine_priority_population], giving the number of people
##'   to vaccinate in each age (row) and priority group (column)
##'
##' @param lag_groups Row indices, corresponding to age
##'   groups in which a lag should be added to the start time of the dose
##'   schedule returned by [vaccine_schedule], if NULL then no lag is added.
##'   Ignored if `lag_groups` is NULL.
##'
##' @param lag_days If `lag_groups` is not NULL then specifies the number of
##'  days to add the start of the dose schedule for the given groups. Ignored
##'  if `lag_groups` is NULL.
##'
##' @param booster_daily_doses_value A vector of booster doses per day.
##'
##' @param booster_groups Which groups of the `priority_population` should be
##'  given a booster, default is all groups; ignored if
##'  `booster_daily_doses_value` is NULL. Sets all rows not corresponding
##'  to given indices to 0.
##'
##' @export
vaccine_schedule_future <- function(start,
                                    daily_doses_value,
                                    mean_days_between_doses,
                                    priority_population,
                                    lag_groups = NULL,
                                    lag_days = NULL,
                                    booster_daily_doses_value = NULL,
                                    booster_groups = 1:19) {

  has_booster <- !is.null(booster_daily_doses_value)

  n_groups <- nrow(priority_population)
  n_priority_groups <- ncol(priority_population)

  if (has_booster) {
    n_doses <- 3L
    n_days <- max(length(daily_doses_value), length(booster_daily_doses_value))
  } else {
    n_doses <- 2L
    n_days <- length(daily_doses_value)
  }

  population_to_vaccinate_mat <-
    array(0, c(n_groups, n_priority_groups, n_doses, n_days))

  population_left <- array(rep(c(priority_population), n_doses),
                           c(n_groups, n_priority_groups, n_doses))

  if (inherits(start, "vaccine_schedule")) {
    for (dose in seq_len(min(n_doses, ncol(start$doses)))) {
      n <- rowSums(start$doses[, dose, ])
      for (i in seq_len(n_priority_groups)) {
        m <- pmin(n, population_left[, i, dose])
        n <- n - m
        population_left[, i, dose] <- population_left[, i, dose] - m
      }
    }

    daily_doses_prev <- apply(start$doses, c(2, 3), sum)
    n_prev <- ncol(daily_doses_prev)
    daily_doses_date <- start$date
  } else {
    daily_doses_prev <- matrix(0, n_doses, 0)
    n_prev <- 0L
    daily_doses_date <- start
  }

  if (nrow(daily_doses_prev) < n_doses) {
    daily_doses_prev <-
      rbind(daily_doses_prev,
            matrix(0, n_doses - nrow(daily_doses_prev),
                   ncol(daily_doses_prev)))
  }
  daily_doses_tt <- cbind(daily_doses_prev, matrix(0, n_doses, n_days))

  population_to_vaccinate_mat <- vaccination_schedule_exec(
    daily_doses_tt, daily_doses_value, population_left,
    population_to_vaccinate_mat, mean_days_between_doses, n_prev, 1:2)
  if (has_booster) {
    ## crude method to zero out any groups not given booster
    booster_population_left <- population_left
    booster_population_left[-booster_groups, , ] <- 0

    population_to_vaccinate_mat <- vaccination_schedule_exec(
      daily_doses_tt, booster_daily_doses_value, booster_population_left,
      population_to_vaccinate_mat, Inf, n_prev, 3)
  }

  doses <- apply(population_to_vaccinate_mat, c(1, 3, 4), sum)

  if (inherits(start, "vaccine_schedule")) {
    if (ncol(start$doses) < n_doses) {
      sdoses <- abind2(start$doses,
                             array(0, c(nrow(start$doses),
                                        n_doses - ncol(start$doses),
                                        nlayer(start$doses))))
    } else {
      sdoses <- start$doses
    }
    doses <- mcstate::array_bind(sdoses, doses)
    dimnames(doses) <- NULL
  }

  schedule <- vaccine_schedule(daily_doses_date, doses, n_doses)

  if (!is.null(lag_groups) || !is.null(lag_days)) {
    if (has_booster) {
      stop("Someone should think about how boost and lag interact")
    }
    if (is.null(lag_days) || is.null(lag_groups)) {
      stop("'lag_days' must be non-NULL iff 'lag_groups' is non_NULL")
    }

    ## we could do this in apply but loop is more readable
    for (i in lag_groups) {
      ## original schedule for group i
      old_schedule_group_i <- schedule$doses[i, , ]
      nr <- nrow(old_schedule_group_i)
      nc <- ncol(old_schedule_group_i)

      ## get start of dose schedule for the group
      start <- which(old_schedule_group_i[1, ] > 0)[1]
      ## catch case when groups are ineligible

      if (!is.na(start)) {
        ## initialize 0 array with same dimensions
        new_schedule_group_i <- array(0, dim(old_schedule_group_i))
        ## add dose schedule to new array after adding the given lag (also
        ##  truncate end by lag amount)
        new_schedule_group_i[, seq.int(start + lag_days, nc)] <-
          old_schedule_group_i[, seq.int(start, nc - lag_days)]
        ## save new schedule
        schedule$doses[i, , ] <- new_schedule_group_i
      }
    }
  }

  schedule
}


##' Create a vaccine schedule for use with [carehomes_parameters]
##'
##' @title Create vaccine schedule
##'
##' @param date A single date, representing the first day that
##'   vaccines will be given
##'
##' @param doses A 3d array of doses representing (1) the model group
##'   (19 rows for the carehomes model), (2) the dose (must be length
##'   2 at present) and (3) time (can be anything nonzero). The values
##'   represent the number of vaccine doses in that group for that
##'   dose for that day. So for `doses[i, j, k]` then it is for the
##'   ith group, the number of jth doses on day `(k - 1) + date`
##'
##' @param n_doses The number of doses in the schedule. Typically (and
##'   by default) this will be 2, but if using booster doses 3 (and in
##'   future we may extend further).
##'
##' @return A `vaccine_schedule` object
##' @export
vaccine_schedule <- function(date, doses, n_doses = 2L) {
  assert_sircovid_date(date)
  assert_scalar(date)

  n_groups <- carehomes_n_groups()

  if (length(dim(doses)) != 3L) {
    stop("Expected a 3d array for 'doses'")
  }
  if (nrow(doses) != n_groups) {
    stop(sprintf("'doses' must have %d rows", n_groups))
  }
  if (ncol(doses) != n_doses) {
    stop(sprintf("'doses' must have %d columns", n_doses))
  }
  if (dim(doses)[[3]] == 0) {
    stop("'doses' must have at least one element in the 3rd dimension")
  }
  if (any(is.na(doses))) {
    stop("'doses' must all be non-NA")
  }
  if (any(doses < 0)) {
    stop("'doses' must all be non-negative")
  }

  ret <- list(date = date, doses = doses, n_doses = n_doses)
  class(ret) <- "vaccine_schedule"
  ret
}


##' Create a historical vaccine schedule from data
##'
##' @title Create historical vaccine schedule
##'
##' @param data A data.frame with columns `date`, `age_band_min`,
##'   `dose1` and `dose2`.
##'
##' @param n_carehomes A vector of length 2 with the number of
##'   carehome workers and residents.
##'
##' @return A [vaccine_schedule] object
##' @export
vaccine_schedule_from_data <- function(data, n_carehomes) {
  assert_is(data, "data.frame")
  required <- c("age_band_min", "date", "dose1", "dose2")
  msg <- setdiff(required, names(data))
  if (length(msg) > 0) {
    stop("Required columns missing from 'data': ",
         paste(squote(msg), collapse = ", "))
  }
  if (length(n_carehomes) != 2) {
    stop("Expected a vector of length 2 for n_carehomes")
  }
  err <- is.na(data$age_band_min) | data$age_band_min %% 5 != 0
  if (any(err)) {
    stop("Invalid values for data$age_band_min: ",
         paste(unique(data$age_band_min[err]), collapse = ", "))
  }
  ## TODO: tidy up later:
  stopifnot(
    all(!is.na(n_carehomes)),
    all(n_carehomes >= 0),
    !is.na(data$date),
    all(data$dose1 >= 0 | is.na(data$dose1)),
    all(data$dose2 >= 0 | is.na(data$dose2)))

  ## First aggregate all the 80+ into one group
  data$date <- as_sircovid_date(data$date)
  data$age_band_min <- pmin(data$age_band_min, 80)
  data <- stats::aggregate(data[c("dose1", "dose2")],
                           data[c("age_band_min", "date")],
                           sum)

  dates <- seq(min(data$date), max(data$date), by = 1)
  age_start <- sircovid_age_bins()$start

  doses <- lapply(c("dose1", "dose2"), function(i)
    stats::reshape(data[c("date", "age_band_min", i)],
                   direction = "wide", timevar = "date",
                   idvar = "age_band_min"))
  stopifnot(identical(dim(doses[[1]]), dim(doses[[2]])))

  ## TODO: add a test for missing days
  i <- match(age_start, doses[[1]]$age_band_min)
  j <- match(dates, sub("^dose[12]\\.", "", names(doses[[1]])))
  doses <- array(
    unlist(lapply(doses, function(d) unname(as.matrix(d)[i, j]))),
    c(length(age_start), length(dates), 2))
  doses <- aperm(doses, c(1, 3, 2))
  doses[is.na(doses)] <- 0

  doses <- vaccine_schedule_add_carehomes(doses, n_carehomes)
  vaccine_schedule(dates[[1]], doses)
}


##' Helper function to create a vaccination schedule that covers data
##' from the past and projects doses into the future based on the last
##' week of vaccination (based on JCVI order using
##' [vaccine_schedule_future]). This function is subject to change.
##'
##' @title Vaccination schedule using data and future
##'
##' @inheritParams vaccine_schedule_from_data
##' @inheritParams vaccine_schedule_future
##' @inheritParams vaccine_priority_population
##'
##' @param end_date The final day in the future to create a schedule
##'   for. After this date the model will assume 0 vaccine doses given
##'   so an overestimate is probably better than an underestimate.
##'
##' @return A [vaccine_schedule] object
##' @export
vaccine_schedule_data_future <- function(data, region, uptake, end_date,
                                         mean_days_between_doses) {
  priority_population <- vaccine_priority_population(region, uptake)
  n_carehomes <- priority_population[18:19, 1]
  schedule_past <- vaccine_schedule_from_data(data, n_carehomes)
  ## then average out last 7 days and project forward by n days
  ## future; we will probably change this to avoid backfill by taking
  ## the days -14..-8
  i <- utils::tail(seq_len(dim(schedule_past$doses)[[3]]), 7)
  mean_doses_last <- sum(schedule_past$doses[, , i], na.rm = TRUE) / length(i)
  end_date <- as_sircovid_date(end_date)
  n_days_future <- end_date -
    (schedule_past$date + dim(schedule_past$doses)[[3]]) + 1L
  daily_doses_future <- rep(round(mean_doses_last), n_days_future)
  vaccine_schedule_future(schedule_past,
                          daily_doses_future,
                          mean_days_between_doses,
                          priority_population)
}


##' @importFrom stats rmultinom
vaccine_schedule_add_carehomes <- function(doses, n_carehomes) {
  doses <- doses[c(seq_len(nrow(doses)), NA, NA), , ]
  doses[is.na(doses)] <- 0

  if (all(n_carehomes == 0)) {
    return(doses)
  }

  ## Impute CHW / CHR - these will be the first vaccinated people
  ## below and above 65 respectively.
  f <- function(target, i_from, i_to) {
    for (i_dose in 1:2) {
      n_t <- colSums(doses[i_from, i_dose, ])
      n <- cumsum(n_t)
      i <- n >= target
      if (!any(i)) {
        k <- length(i)
      } else {
        j <- which(i)[[1]] # the interval that switches
        k <- j - 1
        n_general <- n[[j]] - target
        doses[i_to, i_dose, j] <- target - sum(n_t[seq_len(k)])
        doses[i_from, i_dose, j] <-
          drop(rmultinom(1, n_general, doses[i_from, i_dose, j]))
      }
      doses[i_to, i_dose, seq_len(k)] <- n_t[seq_len(k)]
      doses[i_from, i_dose, seq_len(k)] <- 0
    }
    doses
  }

  age_start <- sircovid_age_bins()$start
  i_chw_from <- which(age_start < 65)
  i_chr_from <- which(age_start >= 65)
  i_chw_to <- 18
  i_chr_to <- 19
  doses <- f(n_carehomes[[1]], i_chw_from, i_chw_to)
  doses <- f(n_carehomes[[2]], i_chr_from, i_chr_to)
  doses
}


##' Create a vaccination scenario
##'
##' @title High-level vaccine scenario creation
##'
##' @param schedule_past A [vaccine_schedule] object corresponding to
##'   previously carried out vaccination.
##'
##' @param doses_future A named vector of vaccine doses to give in the
##'   future. Names must be in ISO date format.
##'
##' @param boosters_future Optional named vector of booster doses to give in
##'   the future. Names must be in ISO date format.
##'
##' @param boosters_prepend_zero If TRUE (default) and `boosters_future` is
##'   not NULL then adds sets booster doses to zero before the first date in
##'   `boosters_future`. This is in contrast to when it is FALSE and the
##'   previous value in `schedule_past` is replicated until the first date in
##'   boosters_future. Note that this should rarely be FALSE as this will
##'   likely lead to duplicating daily doses that are already replicated in
##'   `doses_future`.
##'
##' @inheritParams vaccine_schedule_future
##' @inheritParams vaccine_schedule_data_future
##'
##' @return A [vaccine_schedule] object
##' @export
vaccine_schedule_scenario <- function(schedule_past, doses_future, end_date,
                                      mean_days_between_doses,
                                      priority_population, lag_groups = NULL,
                                      lag_days = NULL,
                                      boosters_future = NULL,
                                      boosters_prepend_zero = TRUE,
                                      booster_groups = 1:19) {

  assert_is(schedule_past, "vaccine_schedule")

  date_end_past <- schedule_past$date + dim(schedule_past$doses)[[3]] - 1L
  i <- utils::tail(seq_len(dim(schedule_past$doses)[[3]]), 7)
  mean_doses_last <- sum(schedule_past$doses[, , i], na.rm = TRUE) / length(i)

  end_date <- as_sircovid_date(end_date)

  if (length(doses_future) == 0 && length(boosters_future) == 0) {
    if (end_date < date_end_past) {
      stop(sprintf(
        "'end_date' must be at least %s (previous end date) but was %s",
        sircovid_date_as_date(date_end_past),
        sircovid_date_as_date(end_date)))
    }
    date_future <- end_date
    booster_daily_doses_value <- NULL
  } else {
    if (length(doses_future) > 0) {
      tmp <- check_doses_boosters_future(doses_future, end_date,
                                         date_end_past)
      date_future <- tmp$date
      doses_future <- tmp$doses
    } else {
      date_future <- end_date
    }
    if (length(boosters_future) > 0) {
      if (boosters_prepend_zero) {
        boosters_future <- c(0, boosters_future)
        names(boosters_future)[1] <-
          as.character(sircovid_date_as_date(date_end_past))
        boosters_future <- boosters_future[!duplicated(boosters_future)]
      }
      tmp <- check_doses_boosters_future(boosters_future, end_date,
                                         date_end_past)

      booster_daily_doses_value <- c(
        rep(mean_doses_last, tmp$date[[1]] - date_end_past),
        rep(unname(tmp$doses), diff(tmp$date)))
    } else {
      booster_daily_doses_value <- NULL
    }
  }

  daily_doses_value <- c(
    rep(mean_doses_last, date_future[[1]] - date_end_past),
    rep(unname(doses_future), diff(date_future)))

  vaccine_schedule_future(schedule_past,
                          daily_doses_value,
                          mean_days_between_doses,
                          priority_population,
                          lag_groups,
                          lag_days,
                          booster_daily_doses_value,
                          booster_groups)
}


vaccination_schedule_exec <- function(daily_doses_tt, daily_doses_value,
                                      population_left,
                                      population_to_vaccinate_mat,
                                      mean_days_between_doses,
                                      n_prev, dose_index) {
  n_priority_groups <- ncol(population_left)
  i1 <- dose_index[[1]]
  i2 <- if (length(dose_index) == 2) dose_index[[2]] else dose_index[[1]]
  for (t in seq_along(daily_doses_value)) {
    tt <- t + n_prev
    tt_dose_1 <- tt - mean_days_between_doses
    if (tt_dose_1 >= 1) {
      ## If we have promised more 2nd doses than we can deliver, we
      ## move our debt forward in time by one day. If doses fluctuate
      ## this will eventually be paid off.
      if (daily_doses_tt[i1, tt_dose_1] > daily_doses_value[t]) {
        daily_doses_tt[i1, tt_dose_1 + 1] <-
          daily_doses_tt[i1, tt_dose_1 + 1] +
          (daily_doses_tt[i1, tt_dose_1] - daily_doses_value[t])
      }
      daily_doses_tt[i2, tt] <- min(daily_doses_value[t],
                                    daily_doses_tt[i1, tt_dose_1])
      daily_doses_tt[i1, tt] <- daily_doses_value[t] - daily_doses_tt[i2, tt]
    } else {
      ## Only distribute first doses
      daily_doses_tt[i2, tt] <- 0
      daily_doses_tt[i1, tt] <- daily_doses_value[t]
    }
    daily_doses_today <- daily_doses_tt[, tt]

    for (dose in dose_index) {
      eligible <- colSums(population_left[, , dose])
      ## Vaccinate the entire of the top priority groups
      n_full_vacc <- findInterval(daily_doses_today[dose], cumsum(eligible))
      if (n_full_vacc > 0) {
        i_full_vacc <- seq_len(n_full_vacc)
        population_to_vaccinate_mat[, i_full_vacc, dose, t] <-
          population_left[, i_full_vacc, dose]
      }

      ## Then partially vaccinate the next priority group, if possible
      if (n_full_vacc < n_priority_groups) {
        if (n_full_vacc == 0) {
          remaining_eligible <- daily_doses_today[dose]
        } else {
          remaining_eligible <- daily_doses_today[dose] -
            cumsum(eligible)[n_full_vacc]
        }
        i_vacc <- n_full_vacc + 1L

        ## Split remaining doses according to age
        population_to_vaccinate_mat[, i_vacc, dose, t] <-
          round(remaining_eligible * population_left[, i_vacc, dose] /
                  sum(population_left[, i_vacc, dose]))
      }

      population_left[, , dose] <- population_left[, , dose] -
        population_to_vaccinate_mat[, , dose, t]
    }
  }

  population_to_vaccinate_mat
}


check_doses_boosters_future <- function(doses, end, end_past) {
  if (is.null(names(doses))) {
    stop(sprintf("'%s' must be named", deparse(substitute(doses))))
  }
  assert_date_string(names(doses), name = sprintf("names(%s)",
                                                  deparse(substitute(doses))))
  doses_future_date <- sircovid_date(names(doses))
  assert_increasing(doses_future_date,
                    name = sprintf("names(%s)",
                                   deparse(substitute(doses))))

  if (last(doses_future_date) > end) {
    stop(sprintf(
      "'end_date' must be at least %s (last %s date) but was %s",
      last(names(doses)),
      deparse(substitute(doses)),
      sircovid_date_as_date(end)))
  }

  if (doses_future_date[[1]] < end_past) {
    message("Trimming vaccination schedule as overlaps with past")
    i <- max(which(doses_future_date < end_past))
    j <- seq(i, length(doses_future_date))
    doses_future_date <- doses_future_date[j]
    doses <- doses[j]
    doses_future_date[[1]] <- end_past
  }

  stopifnot(all(!is.na(doses)))

  date_future <- c(doses_future_date, end)
  names(doses) <- NULL

  list(date = date_future, doses = doses)
}

##' Applies a modifier to severity and transmission parameters
##'
##' @title Modify severity and transmission of variants
##'
##' @param efficacy,efficacy_strain_2 Vaccine efficacy parameters for strains
##'  1 and 2 respectively. Expects a list
##'  with names rel_susceptibility, rel_p_sympt, rel_p_hosp_if_sympt,
##'  rel_infectivity, rel_p_death. Element columns correspond to vaccine strata
##'  and rows to age groups.
##'
##' @param strain_severity_modifier List of modifiers to be applied to efficacy
##'  variables; length should correspond to number of strains and each element
##'  should be a list with same names as `efficacy`
##'
##' @return Returns a list with same length and names as `efficacy` and where
##'  each element has dimensions n_groups x n_strains x n_vacc_strata
##'
##' @export
modify_severity <- function(efficacy, efficacy_strain_2,
                            strain_severity_modifier) {

  expected <- c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt",
                "rel_infectivity", "rel_p_death")

  if (length(efficacy_strain_2) == 0) {
    stopifnot(setequal(names(efficacy), expected))
    n_strain <- 1
  } else {
    stopifnot(
      setequal(names(efficacy), expected),
      setequal(names(efficacy), names(efficacy_strain_2)),
      setequal(lapply(efficacy, dim), lapply(efficacy_strain_2, dim)))
    n_strain <- 4
  }

  n_vacc_strata <- ncol(efficacy[[1]])
  n_groups <- nrow(efficacy[[1]])

  dim <- c(n_groups, n_strain, n_vacc_strata)
  rel_list <- rep(list(array(rep(NA_integer_), dim = dim)), length(expected))
  names(rel_list) <- expected

  for (rel in names(rel_list)) {
    for (s in seq_len(n_strain)) {
      mod <- if (is.null(strain_severity_modifier)) 1 else
        strain_severity_modifier[[s]][[rel]]
      for (v_s in seq_len(n_vacc_strata)) {
        for (g in seq_len(n_groups)) {
          ## If not multistrain then all use same params, otherwise split
          ##  by strains. Strains: 1 (=1), 2(=2), 3(=1->2), 4(=2->1)
          es <- if (n_strain == 1 || s %in% c(1, 4)) efficacy else
            efficacy_strain_2
          new_prob <- es[[rel]][g, v_s] * mod
          ## FIXME - Temporary fixes for when p > 1
          if (new_prob > 1) {
            ## if difference very small or in a class where vaccine hasn't taken
            ##  effect yet, cap at 1
            if (abs(new_prob - 1) < 1e-10 || v_s <= 2) {
              new_prob <- 1
              ## otherwise this is a problem
            } else if (v_s <= 2) {
              browser()
            }
          }
          rel_list[[rel]][g, s, v_s] <- new_prob
        }
      }
    }
  }

  rel_list
}

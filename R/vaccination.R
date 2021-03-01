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

jcvi_prop_to_vaccinate_old <- function(uptake_by_age, 
                                   prop_hcw_by_age,
                                   prop_very_vulnerable_by_age,
                                   prop_underlying_condition_by_age) {
  ## https://www.gov.uk/government/publications/priority-groups-for-coronavirus-covid-19-vaccination-advice-from-the-jcvi-30-december-2020/joint-committee-on-vaccination-and-immunisation-advice-on-priority-groups-for-covid-19-vaccination-30-december-2020#fnref:3
  ## Assuming independance between job (e.g. HCW) and clinical condition
  ## But assuming one is either counted as "clinically extremely vulnerable" or "with underlying health conditions" but not both
  
  n_age_groups <- length(uptake_by_age)
  n_priority_groups <- 12
  p <- matrix(0, n_age_groups, n_priority_groups)
  
  ## JCVI group 1: residents in a care home for older adults and their carers
  j <- 1
  i <- 18:19
  p[i, j] <- 1
  ## JCVI group 2: all those 80 years of age and over and frontline health and social care workers
  j <- 2
  i <- 17
  p[i, j] <- 1
  p[1:(i-1), j] <- prop_hcw_by_age[1:(i-1)]
  ## JCVI group 3: all those 75 years of age and over
  j <- 3
  i <- 16
  p[i, j] <- (1 - prop_hcw_by_age[i])
  ## JCVI group 4: all those 70 years of age and over and clinically extremely vulnerable individuals
  j <- 4
  i <- 15
  p[i, j] <- (1 - prop_hcw_by_age[i])
  p[1:(i-1), j] <- (1 - prop_hcw_by_age[1:(i-1)]) * prop_very_vulnerable_by_age[1:(i-1)]
  ## JCVI group 5: all those 65 years of age and over
  j <- 5
  i <- 14
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i])
  ## JCVI group 6: all individuals aged 16 years to 64 years with underlying health conditions which put them at higher risk of serious disease and mortality
  j <- 6
  i <- 4:13 # starting at 15 not 16 so will need adjusting proportion for this
  p[i, j] <- (1 - prop_hcw_by_age[i]) * prop_underlying_condition_by_age[i]
  ## JCVI group 7: all those 60 years of age and over
  j <- 7
  i <- 13
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i] - prop_underlying_condition_by_age[i])
  ## JCVI group 8: all those 55 years of age and over
  j <- 8
  i <- 12
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i] - prop_underlying_condition_by_age[i])
  ## JCVI group 9: all those 50 years of age and over
  j <- 9
  i <- 11
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i] - prop_underlying_condition_by_age[i])
  
  ## JCVI second phase
  ## JCVI group 10: all those 40-49 years of age and over
  j <- 10
  i <- 9:10
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i] - prop_underlying_condition_by_age[i])
  ## JCVI group 11: all those 30-39 years of age and over
  j <- 11
  i <- 7:8
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i] - prop_underlying_condition_by_age[i])
  ## JCVI group 11: all those 18-29 years of age and over
  j <- 12
  i <- 1:6
  p[i, j] <- (1 - prop_hcw_by_age[i]) * (1 - prop_very_vulnerable_by_age[i] - prop_underlying_condition_by_age[i])
  
  uptake_by_age_mat <- matrix(rep(uptake_by_age, n_priority_groups), nrow = n_age_groups)
  
  p * uptake_by_age_mat
}

jcvi_prop_to_vaccinate <- function(uptake_by_age, 
                                   prop_hcw_by_age,
                                   prop_very_vulnerable_by_age,
                                   prop_underlying_condition_by_age) {
  
  ## https://www.gov.uk/government/publications/priority-groups-for-coronavirus-covid-19-vaccination-advice-from-the-jcvi-30-december-2020/joint-committee-on-vaccination-and-immunisation-advice-on-priority-groups-for-covid-19-vaccination-30-december-2020#fnref:3
  ## Assuming independance between job (e.g. HCW) and clinical condition
  ## But assuming one is either counted as "clinically extremely vulnerable" or "with underlying health conditions" but not both
  
  ## JCVI group 1: residents in a care home for older adults and their carers
  ## JCVI group 2: all those 80 years of age and over and frontline health and social care workers
  ## JCVI group 3: all those 75 years of age and over
  ## JCVI group 4: all those 70 years of age and over and clinically extremely vulnerable individuals
  ## JCVI group 5: all those 65 years of age and over
  ## JCVI group 6: all individuals aged 16 years to 64 years with underlying health conditions which put them at higher risk of serious disease and mortality
  ## JCVI group 7: all those 60 years of age and over
  ## JCVI group 8: all those 55 years of age and over
  ## JCVI group 9: all those 50 years of age and over
  ## JCVI group 10: all those 40-49 years of age and over
  ## JCVI group 11: all those 30-39 years of age and over
  ## JCVI group 11: all those 18-29 years of age and over
  
  n_age_groups <- length(uptake_by_age)
  n_priority_groups <- 12
  p <- matrix(0, n_age_groups, n_priority_groups)
  
  ## Aged base priority list
  jcvi_priority <- list(
    ## the age groups targeted in each priority group (see comments above)
    18:19, 17, 16, 15, 14, NULL, 13, 12, 11, 9:10, 7:8, 1:6 
  )
  
  ## 1. Start with non aged based priority: 
  ## helper function
  add_prop_to_vacc <- function(j, idx, prop_to_vaccinate, p) {
    p[idx, j] <- prop_to_vaccinate[idx]
    p
  }
  
  ## Group 2 includes frontline health and social care workers
  p <- add_prop_to_vacc(j = 2, 
                        idx = 1:(jcvi_priority[[2]]-1), 
                        prop_to_vaccinate = prop_hcw_by_age, 
                        p)
  
  ## Group 4 includes clinically extremely vulnerable individuals
  p <- add_prop_to_vacc(j = 4, 
                        idx = 1:(jcvi_priority[[4]]-1), 
                        prop_to_vaccinate = (1 - prop_hcw_by_age) *
                          prop_very_vulnerable_by_age, 
                        p)
  
  ## Group 6 includes all individuals aged 16 years to 64 years with underlying 
  ## health conditions which put them at higher risk of serious disease and mortality
  p <- add_prop_to_vacc(j = 6, 
                        idx = 4:13, 
                        prop_to_vaccinate = (1 - prop_hcw_by_age) *
                          prop_underlying_condition_by_age, 
                        p)
  
  ## 2. Add aged base priority 
  for(j in seq_along(jcvi_priority)) {
    if (!is.null(jcvi_priority[[j]])) {
      ## discount those already vaccinated as part of non age based priority
      p[jcvi_priority[[j]], j] <- 1 - rowSums(p)[jcvi_priority[[j]]]
    }
  }
  
  ## 3. Account for uptake
  uptake_by_age_mat <- matrix(rep(uptake_by_age, n_priority_groups), nrow = n_age_groups)
  
  p * uptake_by_age_mat
}


n_to_vaccinate <- function(prop_to_vaccinate, region) {
  pop_by_age <- sircovid:::carehomes_parameters(1, region)$N_tot
  pop_by_age_mat <- matrix(rep(pop_by_age, ncol(prop_to_vaccinate)),
                           nrow = nrow(prop_to_vaccinate))
  prop_to_vaccinate * pop_by_age_mat
}


jcvi_n_to_vaccinate <- function(uptake_by_age, 
                                prop_hcw_by_age,
                                prop_very_vulnerable_by_age,
                                prop_underlying_condition_by_age,
                                region) {
  p <- jcvi_prop_to_vaccinate(uptake_by_age, 
                              prop_hcw_by_age,
                              prop_very_vulnerable_by_age,
                              prop_underlying_condition_by_age)
  n_to_vaccinate(p, region)
}


get_dose_schedule <- function(daily_doses, n, mean_days_between_doses) {
  ## n matrix with dimensions n_age_groups, n_priority_groups
  ## as obtained from jcvi_n_to_vaccinate function
  n_age_groups <- nrow(n)
  n_priority_groups <- ncol(n)
  
  # first doses
  first_daily_doses <- 0 * daily_doses
  n_to_vaccinate_mat <- array(0, dim = c(n_age_groups, n_priority_groups, length(daily_doses)))
  n_to_vaccinate <- array(0, dim = c(n_age_groups, length(daily_doses)))
  n_left <- n
  # second doses
  second_daily_doses <- 0 * daily_doses
  n_to_vaccinate_second_dose <- array(0, dim = c(n_age_groups, length(daily_doses)))
  
  for (t in seq_along(daily_doses)) {
    
    ## split doses between first and second doses
    if (t <= mean_days_between_doses) { # only distribute first doses
      second_daily_doses[t] <- 0  
      first_daily_doses[t] <- daily_doses[t]
    } else { # prioritise second doses
      second_daily_doses[t] <- min(daily_doses[t], first_daily_doses[t - mean_days_between_doses])
      first_daily_doses[t] <- daily_doses[t]  - second_daily_doses[t] 
    }
    
    ## allocate first doses
    eligible <- colSums(n_left)
    ## vaccinate fully the top priority groups
    n_full_vacc <- findInterval(first_daily_doses[t], cumsum(eligible)) # vaccinate everyone in these groups
    if (n_full_vacc > 0) {
      n_to_vaccinate_mat[, 1:n_full_vacc, t] <- n_left[,1:n_full_vacc] 
    }
    ## then go down to next priority group
    if (n_full_vacc < n_priority_groups) {
      if(n_full_vacc == 0) {
        remaining_eligible <- first_daily_doses[t]
      } else {
        remaining_eligible <- first_daily_doses[t] - cumsum(eligible)[n_full_vacc]
      }
      ## split remaining doses according to age
      n_to_vaccinate_mat[, n_full_vacc + 1, t] <- round(remaining_eligible * n_left[,n_full_vacc + 1] / sum(n_left[,n_full_vacc + 1] ))
    }
    
    n_left <- n_left - n_to_vaccinate_mat[, , t]
    
    n_to_vaccinate[, t] <- rowSums(n_to_vaccinate_mat[, , t])
    
    ## allocate second doses 12 weeks later
    if (t > mean_days_between_doses) { # also do second doses
      n_to_vaccinate_second_dose[, t] <- n_to_vaccinate[, t - mean_days_between_doses]
    }
    
  } # end loop over time
  
  list(n_to_vaccinate = n_to_vaccinate,
       n_to_vaccinate_second_dose = n_to_vaccinate_second_dose)
  
}

##' Compute "Rt" for a single simulated trajectory and parameter set.
##'
##' @title Compute "Rt"
##'
##' @param step A vector of steps that the model was run over
##'
##' @param S A 19 x steps matrix of "S" compartment counts
##'
##' @param p A [carehomes_parameters()] object
##'
##' @return A list with elements `step`, `beta`, `eff_Rt_all`,
##'   `eff_Rt_general`, `Rt_all` and `Rt_general`
##'
##' @export
carehomes_Rt <- function(step, S, p) {
  if (nrow(S) != ncol(p$rel_susceptibility) * nrow(p$m)) {
    stop(sprintf(
      "Expected 'S' to have %d rows, following transmission matrix",
      nrow(p$m)))
  }
  if (ncol(S) != length(step)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'step'",
                 length(step)))
  }

  beta <- sircovid_parameters_beta_expand(step, p$beta_step)
  mean_duration <- carehomes_Rt_mean_duration(step, p)
  max_strain_multiplier <- max(p$strain_transmission)

  calculate_ev <-
    function(t, S, beta, mean_duration, max_strain_multiplier, drop_carehomes) {
    ## Next-Generation-Matrix
    m <- p$m
    ages <- seq_len(p$n_age_groups)
    ch <- seq(to = p$n_groups, length.out = 2)
    m[ages, ] <- beta[t] * m[ages, ]
    m[ch, ages] <- beta[t] * m[ch, ages]

    ## when several vaccination groups,
    ## need to take the weighted means of the S
    ## (weights given by rel_susceptibility)
    S_mat <- matrix(S[, t],
                    nrow = p$n_groups, ncol = ncol(p$rel_susceptibility))
    S_weighted <- rowSums(S_mat * p$rel_susceptibility)

    if (dim(mean_duration)[2] > 1) {
      mean_duration_weighted <- apply(mean_duration[, , t], 1, mean)
    } else {
      mean_duration_weighted <- drop(mean_duration[, , t])
    }

    ## In a multistrain model R0 is the max of R0 across strains
    ngm <- outer(mean_duration_weighted, S_weighted) * m * max_strain_multiplier

    ## Care home workers (CHW) and residents (CHR) in last two rows
    ## and columns
    if (drop_carehomes) {
      i <- seq_len(nrow(ngm) - 2L)
      ngm <- ngm[i, i]
    }
    ev <- eigen(ngm)$values
    ev[Im(ev) != 0] <- NA
    ev <- as.numeric(ev)
    max(ev, na.rm = TRUE)
  }

  t <- seq_along(step)
  eff_Rt_all <- vnapply(t, calculate_ev, S,
                        beta = beta,
                        mean_duration = mean_duration,
                        max_strain_multiplier = max_strain_multiplier,
                        drop_carehomes = FALSE)
  eff_Rt_general <- vnapply(t, calculate_ev, S,
                            beta = beta,
                            mean_duration = mean_duration,
                            max_strain_multiplier = max_strain_multiplier,
                            drop_carehomes = TRUE)
  N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
  N_tot_all_vacc_groups <- N_tot_non_vacc
  for (i in seq(2, ncol(p$rel_susceptibility))) {
    N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                   0 * N_tot_non_vacc)
  }
  Rt_all <- vnapply(t, calculate_ev, N_tot_all_vacc_groups,
                    beta = beta,
                    mean_duration = mean_duration,
                    max_strain_multiplier = max_strain_multiplier,
                    drop_carehomes = FALSE)
  Rt_general <- vnapply(t, calculate_ev, N_tot_all_vacc_groups,
                        beta = beta,
                        mean_duration = mean_duration,
                        max_strain_multiplier = max_strain_multiplier,
                        drop_carehomes = TRUE)

  list(step = step,
       date = step * p$dt,
       beta = beta,
       eff_Rt_all = eff_Rt_all,
       eff_Rt_general = eff_Rt_general,
       Rt_all = Rt_all,
       Rt_general = Rt_general)
}

## Here we expect 'S' in order:
##
##   state x sample x step
##
## We expect 'pars' to be a list along sample (or a shared parameter set)
## We expect 'step' to be a vector along step

##' Compute "Rt" for a set of simulated trajectories (e.g., the result
##' of [dust::dust_iterate()], [mcstate::pmcmc()] or
##' [mcstate::pmcmc_predict()]. The trajectories may or may not share
##' parameters.
##'
##' @title Compute Rt for a set of trajectories
##'
##' @param step A vector of steps
##'
##' @param S A 3d (19 x n trajectories x n steps) array of "S"
##'   compartment counts
##'
##' @param pars Either a single [carehomes_parameters()] object
##'   (shared parameters) or an unnamed list of
##'   [carehomes_parameters()] objects, the same length as `ncol(S)`.
##'
##' @param initial_step_from_parameters If `TRUE`, then `step[[1]]` is
##'   replaced by the value of `initial_step` from the parameters.
##'   This is usually what you want.
##'
##' @param shared_parameters Should `pars` be treated as a single
##'   shared list? Leave as `NULL` to detect automatically, set to
##'   `TRUE` or `FALSE` to force it to be interpreted one way or the
##'   other which may give more easily interpretable error messages.
##'
##' @return As for [carehomes_Rt()], except that every element is a
##'   matrix, not a vector.
##'
##' @export
carehomes_Rt_trajectories <- function(step, S, pars,
                                      initial_step_from_parameters = TRUE,
                                      shared_parameters = NULL) {
  calculate_Rt_trajectories(carehomes_Rt, step, S, pars,
                            initial_step_from_parameters, shared_parameters)
}


carehomes_Rt_mean_duration <- function(step, pars) {
  dt <- pars$dt

  matricise <- function(vect, n_col) {
    matrix(rep(vect, n_col), ncol = n_col, byrow = FALSE)
  }

  n_vacc_classes <- ncol(pars$rel_susceptibility)

  n_groups <- pars$n_groups

  n_time_steps <-
    length(sircovid_parameters_beta_expand(step, pars$p_H_step))

  ## TODO: This is not correct for the initial transition from
  ## vaccination
  vacc_prog_before_infectious <- function(group_i) {
    vacc_prog_rate <- pars$vaccine_progression_rate_base[group_i, ]

    ## Q[i, j] gives probability of progression from vaccine stage i to
    ## vaccine stage j in one time step
    Q <- diag(exp(-vacc_prog_rate * dt))
    for (i in seq_len(n_vacc_classes - 1)) {
      Q[i, i + 1] <- 1 - Q[i, i]
    }
    Q[n_vacc_classes, 1] <- 1 - Q[n_vacc_classes, n_vacc_classes]

    ## probability of E progression in one time step
    p_EE <- 1 - exp(-pars$gamma_E * dt)

    ## A[i, j] gives the probability that an individual who begins an E stage in
    ## vaccine stage i, exits that E stage in vaccine stage j. Note that A =
    ## (sum_{k = 0}^Inf ((1 - p_EE) * Q) ^ k) * p_EE * Q. Also note that for a
    ## square matrix B, sum_{k = 0}^Inf B^k = (I - B)^-1.
    A <- (solve(diag(n_vacc_classes) - (1 - p_EE) * Q) %*% (p_EE * Q))

    ## Note we need to account for there being s_E stages in E, and also that
    ## individuals can have a vaccine progression in the same step that they get
    ## infected (hence Q appearing below).
    out <- Q %*% matrix_pow(A, pars$s_E)
    out
  }

  if (n_vacc_classes > 1) {
    ## V[i, j, k] gives the probability that an individual in group k who is
    ## infected when in vaccine stage i exits the E class in vaccine stage j
    V <- array(sapply(seq_len(pars$n_groups), vacc_prog_before_infectious),
             dim = c(n_vacc_classes, n_vacc_classes, pars$n_groups))
  }

  mat_multi_by_group <- function(p, V) {
    out <- sapply(seq_len(pars$n_groups),
                  function(i) {
                    V[, , i] %*% p[i, ]})
    out
  }

  p_C <- matricise(pars$p_C, n_vacc_classes)
  p_C <- p_C * pars$rel_p_C
  if (n_vacc_classes > 1) {
    p_C <- t(mat_multi_by_group(p_C, V))
  }
  p_C <- outer(p_C, rep(1, n_time_steps))

  p_H <- matricise(pars$psi_H, n_vacc_classes) *
    pars$rel_p_H
  if (n_vacc_classes > 1) {
    p_H <- t(mat_multi_by_group(p_H, V))
  }
  p_H <- outer(p_H,
    sircovid_parameters_beta_expand(step, pars$p_H_step))

  p_ICU_hosp <- outer(matricise(pars$psi_ICU_hosp, n_vacc_classes),
                    sircovid_parameters_beta_expand(step, pars$p_ICU_hosp_step))
  p_death_ICU <- outer(matricise(pars$psi_death_ICU, n_vacc_classes),
                   sircovid_parameters_beta_expand(step, pars$p_death_ICU_step))
  p_death_hosp_D <- outer(matricise(pars$psi_death_hosp_D, n_vacc_classes),
                sircovid_parameters_beta_expand(step, pars$p_death_hosp_D_step))
  p_death_stepdown <- outer(matricise(pars$psi_death_stepdown, n_vacc_classes),
              sircovid_parameters_beta_expand(step, pars$p_death_stepdown_step))
  p_death_comm <- outer(matricise(pars$psi_death_comm, n_vacc_classes),
                  sircovid_parameters_beta_expand(step, pars$p_death_comm_step))

  p_hosp_R <- p_C * p_H * (1 - p_death_comm) *
    (1 - p_ICU_hosp) * (1 - p_death_hosp_D)
  p_hosp_D <- p_C * p_H * (1 - p_death_comm) *
    (1 - p_ICU_hosp) * p_death_hosp_D
  p_ICU_S_R <- p_C * p_H * (1 - p_death_comm) *
    p_ICU_hosp * (1 - p_death_ICU) * (1 - p_death_stepdown)
  p_ICU_S_D <- p_C * p_H * (1 - p_death_comm) *
    p_ICU_hosp * (1 - p_death_ICU) * p_death_stepdown
  p_ICU_D <- p_C * p_H * (1 - p_death_comm) *
    p_ICU_hosp * p_death_ICU

  ## TODO: would be nice if it's possibly to name these subcomponents
  ## to make the calculation clearer.
  mean_duration <- (1 - p_C) * pars$s_A /
    (1 - exp(- dt * pars$gamma_A)) +
    p_C * pars$s_C / (1 - exp(- dt * pars$gamma_C))

  mean_duration <- mean_duration +
    pars$comm_D_transmission * p_C * p_H *
    p_death_comm * pars$s_G_D / (1 - exp(- dt * pars$gamma_G_D))

  mean_duration <- mean_duration +
    pars$hosp_transmission * (
      p_hosp_R * pars$s_hosp_R / (1 - exp(- dt * pars$gamma_hosp_R)) +
      p_hosp_D * pars$s_hosp_D / (1 - exp(- dt * pars$gamma_hosp_D)) +
      (p_ICU_S_R + p_ICU_S_D + p_ICU_D) * pars$s_triage /
      (1 - exp(- dt * pars$gamma_triage))) +
    pars$ICU_transmission * (
      p_ICU_S_R * pars$s_ICU_S_R / (1 - exp(- dt * pars$gamma_ICU_S_R)) +
      p_ICU_S_D * pars$s_ICU_S_D / (1 - exp(- dt * pars$gamma_ICU_S_D)) +
      p_ICU_D * pars$s_ICU_D / (1 - exp(- dt * pars$gamma_ICU_D)))

  dt * mean_duration
}


## This part holds over all possible Rt calculations, so I've factored
## it out here; when we implement this for the basic model this will
## remain unchanged.  However, I am leaving it in this
## carehomes-specific file until we do add a new model or port it.
calculate_Rt_trajectories <- function(calculate_Rt, step, S, pars,
                                      initial_step_from_parameters,
                                      shared_parameters) {
  if (length(dim(S)) != 3) {
    stop("Expected a 3d array of 'S'")
  }

  shared_parameters <- shared_parameters %||% !is.null(names(pars))
  if (shared_parameters) {
    if (is.null(names(pars))) {
      stop("If using shared parameters, expected a named list for 'pars'")
    }
    pars <- rep(list(pars), ncol(S))
  } else {
    if (!is.null(names(pars))) {
      stop("If not using shared parameters, expected a unnamed list for 'pars'")
    }
    if (length(pars) != ncol(S)) {
      stop(sprintf(
        "Expected 2nd dimension of 'S' to have length %d, following 'pars'",
        length(pars)))
    }
  }

  if (dim(S)[[3]] != length(step)) {
    stop(sprintf(
      "Expected 3rd dimension of 'S' to have length %d, following 'step'",
      length(step)))
  }

  calculate_rt_one_trajectory <- function(i) {
    if (initial_step_from_parameters) {
      step[[1L]] <- pars[[i]]$initial_step
    }
    calculate_Rt(step, S[, i, ], pars[[i]])
  }

  res <- lapply(seq_along(pars), calculate_rt_one_trajectory)

  ## These are stored in a list-of-lists and we convert to a
  ## list-of-matrices here
  collect <- function(nm) {
    matrix(unlist(lapply(res, "[[", nm)), length(step), length(res))
  }
  nms <- names(res[[1]])
  ret <- set_names(lapply(nms, collect), nms)
  class(ret) <- "Rt_trajectories"
  ret
}

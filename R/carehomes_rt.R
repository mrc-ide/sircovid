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

  calculate_ev <- function(t, S, drop_carehomes) {
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

    ngm <- outer(mean_duration[, t], S_weighted) * m

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
  eff_Rt_all <- vnapply(t, calculate_ev, S, drop_carehomes = FALSE)
  eff_Rt_general <- vnapply(t, calculate_ev, S, drop_carehomes = TRUE)
  N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
  N_tot_all_vacc_groups <- N_tot_non_vacc
  for (i in seq(2, ncol(p$rel_susceptibility))) {
    N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                   0 * N_tot_non_vacc)
  }
  Rt_all <- vnapply(t, calculate_ev, N_tot_all_vacc_groups,
                    drop_carehomes = FALSE)
  Rt_general <- vnapply(t, calculate_ev, N_tot_all_vacc_groups,
                        drop_carehomes = TRUE)

  list(step = step,
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

  p_asympt <- pars$p_asympt
  p_sympt_ILI <- pars$p_sympt_ILI
  p_hosp_ILI <- outer(pars$psi_hosp_ILI,
                    sircovid_parameters_beta_expand(step, pars$p_hosp_ILI_step))
  p_ICU_hosp <- outer(pars$psi_ICU_hosp,
                    sircovid_parameters_beta_expand(step, pars$p_ICU_hosp_step))
  p_death_ICU <- outer(pars$psi_death_ICU,
                   sircovid_parameters_beta_expand(step, pars$p_death_ICU_step))
  p_death_hosp_D <- outer(pars$psi_death_hosp_D,
                sircovid_parameters_beta_expand(step, pars$p_death_hosp_D_step))
  p_death_stepdown <- outer(pars$psi_death_stepdown,
              sircovid_parameters_beta_expand(step, pars$p_death_stepdown_step))
  p_death_comm <- outer(pars$psi_death_comm,
                  sircovid_parameters_beta_expand(step, pars$p_death_comm_step))

  p_mild <- (1 - p_asympt) * (1 - p_sympt_ILI)
  p_ILI <- (1 - p_asympt) * p_sympt_ILI

  p_hosp_R <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    (1 - p_ICU_hosp) * (1 - p_death_hosp_D)
  p_hosp_D <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    (1 - p_ICU_hosp) * p_death_hosp_D
  p_ICU_S_R <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    p_ICU_hosp * (1 - p_death_ICU) * (1 - p_death_stepdown)
  p_ICU_S_D <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    p_ICU_hosp * (1 - p_death_ICU) * p_death_stepdown
  p_ICU_D <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    p_ICU_hosp * p_death_ICU

  ## TODO: would be nice if it's possibly to name these subcomponents
  ## to make the calculation clearer.
  mean_duration <- p_asympt * pars$s_asympt /
    (1 - exp(- dt * pars$gamma_asympt)) +
    p_mild * pars$s_mild / (1 - exp(- dt * pars$gamma_mild)) +
    p_ILI * pars$s_ILI / (1 - exp(- dt * pars$gamma_ILI))

  mean_duration <- mean_duration +
    pars$comm_D_transmission * p_ILI * p_hosp_ILI *
    p_death_comm * pars$s_comm_D / (1 - exp(- dt * pars$gamma_comm_D))

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
  set_names(lapply(nms, collect), nms)
}

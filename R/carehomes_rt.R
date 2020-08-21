## Compute Rt for one trajectory, given the number of susceptibles
## (S), the beta over time (same length as S) and other carehomes
## parameters.  This will be the result of one simulation.
carehomes_Rt <- function(step, S, p) {
  if (nrow(S) != nrow(p$m)) {
    stop(sprintf(
      "Expected 'S' to have %d rows, following transmission matrix",
      nrow(p$m)))
  }
  if (ncol(S) != length(step)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'step'",
                 length(step)))
  }

  ## Expand out beta following the model
  beta <- p$beta_step[pmin(step, length(p$beta_step))]
  mean_duration <- carehomes_Rt_mean_duration(p)

  calculate_ev <- function(t, drop_carehomes) {
    ## Next-Generation-Matrix
    ngm <- beta[t] * outer(mean_duration, S[, t]) * p$m

    ## Care home workers (CHW) and residents (CHR) in last two
    ## rows/column
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
  eff_Rt_all <- vnapply(t, calculate_ev, drop_carehomes = FALSE)
  eff_Rt_general <- vnapply(t, calculate_ev, drop_carehomes = TRUE)

  list(step = step,
       beta = beta,
       eff_Rt_all = eff_Rt_all,
       eff_Rt_general = eff_Rt_general,
       Rt_all = eff_Rt_all[[1]] / beta[[1]] * beta,
       Rt_general = eff_Rt_general[[1]] / beta[[1]] * beta)
}


## Here we expect 'S' in order:
##
##   state / sample / step
##
## We expect 'pars' to be a list along sample (or a shared parameter set)
## We expect 'step' to be a vector along step
carehomes_Rt_trajectories <- function(step, S, pars,
                                      initial_step_from_parameters = TRUE,
                                      shared_parameters = NULL) {
  calculate_Rt_trajectories(carehomes_Rt, step, S, pars,
                            initial_step_from_parameters, shared_parameters)
}


carehomes_Rt_mean_duration <- function(pars) {
  dt <- pars$dt
  p_ICU_hosp  <- pars$p_ICU_hosp
  p_ICU_hosp <- pars$p_ICU_hosp
  p_asympt <- pars$p_asympt
  p_death_ICU <- pars$p_death_ICU
  p_death_comm <- pars$p_death_comm
  p_death_hosp_D <- pars$p_death_hosp_D
  p_hosp_ILI <- pars$p_hosp_ILI
  p_sympt_ILI <- pars$p_sympt_ILI

  p_mild <- (1 - p_asympt) * (1 - p_sympt_ILI)
  p_ILI <- (1 - p_asympt) * p_sympt_ILI

  p_hosp_R <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    (1 - p_ICU_hosp) * (1 - p_death_hosp_D)
  p_hosp_D <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    (1 - p_ICU_hosp) * p_death_hosp_D
  p_ICU_R <- p_ILI * p_hosp_ILI * (1 - p_death_comm) *
    p_ICU_hosp * (1 - p_death_ICU)
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
      p_hosp_R * pars$s_hosp_R/(1 - exp(- dt * pars$gamma_hosp_R)) +
      p_hosp_D * pars$s_hosp_D/(1 - exp(- dt * pars$gamma_hosp_D)) +
      (p_ICU_R + p_ICU_D) * pars$s_triage /
      (1 - exp(- dt * pars$gamma_triage))) +
    pars$ICU_transmission * (
      p_ICU_R * pars$s_ICU_R / (1 - exp(- dt * pars$gamma_ICU_R)) +
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

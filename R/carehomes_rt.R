##' Compute "Rt" for a single simulated trajectory and parameter set.
##'
##' @title Compute "Rt"
##'
##' @param step A vector of steps that the model was run over
##'
##' @param S A (n groups x n vaccine classes) x steps matrix of "S"
##'   compartment counts
##'
##' @param p A [carehomes_parameters()] object
##'
##' @param prob_strain A (n groups x n strains) x n steps matrix of
##'   "prob_strain" outputs from the model. For a 2 strain model for example,
##'   `prob_strain[1, j]` and `prob_strain[n_groups + 1, j]` should give, for the
##'   j^th time step, the probabilities that new infections
##'   in group 1 are of strains 1 and 2 respectively.
##'   The default is `NULL`, but it must
##'   be specified if there is more than one strain
##'
##' @return A list with elements `step`, `beta`, `eff_Rt_all`,
##'   `eff_Rt_general`, `Rt_all` and `Rt_general`
##'
##' @export
carehomes_Rt <- function(step, S, p, prob_strain = NULL) {
  if (nrow(S) != ncol(p$rel_susceptibility) * nrow(p$m)) {
    stop(sprintf(
      "Expected 'S' to have %d rows = %d groups x %d vaccine classes",
      p$n_groups * ncol(p$rel_susceptibility),
      p$n_groups,
      ncol(p$rel_susceptibility)))
  }
  if (ncol(S) != length(step)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'step'",
                 length(step)))
  }
  if (is.null(prob_strain)) {
    if (length(p$strain_transmission) > 1) {
      stop("Expected prob_strain input because there is more than one strain")
    } else {
      prob_strain <- array(1, c(p$n_groups, length(step)))
    }
  } else {
    if (nrow(prob_strain) != length(p$strain_transmission) * p$n_groups) {
      stop(sprintf(
        "Expected 'prob_strain' to have %d rows = %d groups x %d strains",
        p$n_groups * length(p$strain_transmission),
        p$n_groups,
        length(p$strain_transmission)))
    }
    if (ncol(prob_strain) != length(step)) {
      stop(sprintf(
        "Expected 'prob_strain' to have %d columns, following 'step'",
                   length(step)))
    }
  }

  ### here mean_duration accounts for relative infectivity of
  ### different infection / vaccination stages
  beta <- sircovid_parameters_beta_expand(step, p$beta_step)
  mean_duration <- carehomes_Rt_mean_duration_weighted_by_infectivity(step, p)

  n_vacc_classes <- ncol(p$rel_susceptibility)

  calculate_ev <-
    function(t, S, prob_strain, beta, mean_duration, drop_carehomes) {
    ## Next-Generation-Matrix
    m <- p$m
    ages <- seq_len(p$n_age_groups)
    ch <- seq(to = p$n_groups, length.out = 2)
    m[ages, ] <- beta[t] * m[ages, ]
    m[ch, ages] <- beta[t] * m[ch, ages]

    m_extended <- matrix(t(matrix(m, p$n_groups, p$n_groups * n_vacc_classes)),
                         p$n_groups * n_vacc_classes,
                         p$n_groups * n_vacc_classes,
                         byrow = TRUE)

    S_weighted <- S[, t] * c(p$rel_susceptibility)

    prob_strain_mat <- matrix(prob_strain[, t],
                              nrow = p$n_groups,
                              ncol = length(p$strain_transmission))
    weighted_strain_multiplier <- prob_strain_mat %*% p$strain_transmission

    ngm <- outer(c(mean_duration[, , t] *
                     array(weighted_strain_multiplier,
                           c(p$n_groups, n_vacc_classes))),
                 S_weighted) * m_extended

    ## Care home workers (CHW) and residents (CHR) in last two rows
    ## and columns, remove for each vaccine class
    if (drop_carehomes) {
      i_CHR <- seq(p$n_groups, nrow(ngm), by = p$n_groups)
      i_CHW <- i_CHR - 1
      i_gen <- seq_len(nrow(ngm))[-c(i_CHW, i_CHR)]
      ngm <- ngm[i_gen, i_gen]
    }
    ev <- eigen(ngm)$values
    ev[Im(ev) != 0] <- NA
    ev <- as.numeric(ev)
    max(ev, na.rm = TRUE)
  }

  t <- seq_along(step)
  eff_Rt_all <- vnapply(t, calculate_ev, S,
                        prob_strain = prob_strain,
                        beta = beta,
                        mean_duration = mean_duration,
                        drop_carehomes = FALSE)
  eff_Rt_general <- vnapply(t, calculate_ev, S,
                            prob_strain = prob_strain,
                            beta = beta,
                            mean_duration = mean_duration,
                            drop_carehomes = TRUE)
  N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
  N_tot_all_vacc_groups <- N_tot_non_vacc
  if (n_vacc_classes > 1) {
    for (i in 2:n_vacc_classes) {
      N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                     0 * N_tot_non_vacc)
    }
  }
  Rt_all <- vnapply(t, calculate_ev, N_tot_all_vacc_groups,
                    prob_strain = prob_strain,
                    beta = beta,
                    mean_duration = mean_duration,
                    drop_carehomes = FALSE)
  Rt_general <- vnapply(t, calculate_ev, N_tot_all_vacc_groups,
                        prob_strain = prob_strain,
                        beta = beta,
                        mean_duration = mean_duration,
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
##' @param S A 3d ((n groups x n vaccine classes) x n trajectories x n steps)
##'   array of "S" compartment counts
##'
##' @param pars Either a single [carehomes_parameters()] object
##'   (shared parameters) or an unnamed list of
##'   [carehomes_parameters()] objects, the same length as `ncol(S)`.
##'
##' @param prob_strain A 3d ((n groups x n strains) x n trajectories x n steps)
##'   array of "prob_strain" model outputs. Default is `NULL`, but it must be
##'   specified if there is more than one strain.
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
carehomes_Rt_trajectories <- function(step, S, pars, prob_strain = NULL,
                                      initial_step_from_parameters = TRUE,
                                      shared_parameters = NULL) {
  calculate_Rt_trajectories(carehomes_Rt, step, S, pars, prob_strain,
                            initial_step_from_parameters, shared_parameters)
}


carehomes_Rt_mean_duration_weighted_by_infectivity <- function(step, pars) {
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

    ## Note we need to account for there being k_E stages in E, and also that
    ## individuals can have a vaccine progression in the same step that they get
    ## infected (hence Q appearing below).
    out <- Q %*% matrix_pow(A, pars$k_E)
    out
  }

  if (n_vacc_classes > 1) {
    ## V[i, j, k] gives the probability that an individual in group k who is
    ## infected when in vaccine stage i exits the E class in vaccine stage j
    V <- vapply(seq_len(pars$n_groups), vacc_prog_before_infectious,
                array(0, c(n_vacc_classes, n_vacc_classes)))
  }

  ## compute probabilities of different pathways

  p_C <- matricise(pars$p_C, n_vacc_classes) * pars$rel_p_sympt
  p_C <- outer(p_C, rep(1, n_time_steps))

  p_H <- matricise(pars$psi_H, n_vacc_classes) * pars$rel_p_hosp_if_sympt
  p_H <- outer(p_H, sircovid_parameters_beta_expand(step, pars$p_H_step))

  p_ICU <- outer(matricise(pars$psi_ICU, n_vacc_classes),
                 sircovid_parameters_beta_expand(step, pars$p_ICU_step))
  p_ICU_D <- outer(matricise(pars$psi_ICU_D, n_vacc_classes),
                   sircovid_parameters_beta_expand(step, pars$p_ICU_D_step))
  p_H_D <- outer(matricise(pars$psi_H_D, n_vacc_classes),
                sircovid_parameters_beta_expand(step, pars$p_H_D_step))
  p_W_D <- outer(matricise(pars$psi_W_D, n_vacc_classes),
                 sircovid_parameters_beta_expand(step, pars$p_W_D_step))
  p_G_D <- outer(matricise(pars$psi_G_D, n_vacc_classes),
                 sircovid_parameters_beta_expand(step, pars$p_G_D_step))

  prob_H_R <- p_C * p_H * (1 - p_G_D) *
    (1 - p_ICU) * (1 - p_H_D)
  prob_H_D <- p_C * p_H * (1 - p_G_D) *
    (1 - p_ICU) * p_H_D
  prob_ICU_W_R <- p_C * p_H * (1 - p_G_D) *
    p_ICU * (1 - p_ICU_D) * (1 - p_W_D)
  prob_ICU_W_D <- p_C * p_H * (1 - p_G_D) *
    p_ICU * (1 - p_ICU_D) * p_W_D
  prob_ICU_D <- p_C * p_H * (1 - p_G_D) *
    p_ICU * p_ICU_D

  ## Compute mean duration (in time steps) of each stage of infection,
  ## weighed by probability of going through that stage
  ## and by relative infectivity of that stage

  ## Note the mean duration (in time steps) of a compartment for
  ## a discretised Erlang(k, gamma) is k / (1 - exp(dt * gamma))

  mean_duration_I_A <- (1 - p_C) * pars$k_A / (1 - exp(- dt * pars$gamma_A))

  mean_duration_I_C <- p_C * pars$k_C / (1 - exp(- dt * pars$gamma_C))

  mean_duration_G_D <- pars$G_D_transmission * p_C * p_H *
    p_G_D * pars$k_G_D / (1 - exp(- dt * pars$gamma_G_D))

  mean_duration_hosp <- pars$hosp_transmission * (
    prob_H_R * pars$k_H_R / (1 - exp(- dt * pars$gamma_H_R)) +
      prob_H_D * pars$k_H_D / (1 - exp(- dt * pars$gamma_H_D)) +
      (prob_ICU_W_R + prob_ICU_W_D + prob_ICU_D) * pars$k_ICU_pre /
      (1 - exp(- dt * pars$gamma_ICU_pre)))

  mean_duration_icu <- pars$ICU_transmission * (
    prob_ICU_W_R * pars$k_ICU_W_R / (1 - exp(- dt * pars$gamma_ICU_W_R)) +
      prob_ICU_W_D * pars$k_ICU_W_D / (1 - exp(- dt * pars$gamma_ICU_W_D)) +
      prob_ICU_D * pars$k_ICU_D / (1 - exp(- dt * pars$gamma_ICU_D)))

  mean_duration <- mean_duration_I_A + mean_duration_I_C + mean_duration_G_D +
    mean_duration_hosp + mean_duration_icu

  ## Account for different infectivity levels depending on vaccination stage

  mean_duration <- mean_duration *
    outer(pars$rel_infectivity, rep(1, n_time_steps))

  ## Multiply by dt to convert from time steps to days
  mean_duration <- dt * mean_duration

  ## mean_duration[i, j, k] represents mean duration at step k of age group i
  ## leaving the E compartment in vaccine stage j, we need to output for leaving
  ## the S compartment in vaccine stage j, so we calculate this here
  if (n_vacc_classes > 1) {
    out <- array(0, dim(mean_duration))
    for (i in seq_len(pars$n_groups)) {
      for (j in seq_len(n_time_steps)) {
        out[i, , j] <- V[, , i] %*% mean_duration[i, , j]
      }
    }
  } else {
    out <- mean_duration
  }

  out
}


## This part holds over all possible Rt calculations, so I've factored
## it out here; when we implement this for the basic model this will
## remain unchanged.  However, I am leaving it in this
## carehomes-specific file until we do add a new model or port it.
calculate_Rt_trajectories <- function(calculate_Rt, step, S, pars, prob_strain,
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

  if (!is.null(prob_strain)) {
    if (length(dim(prob_strain)) != 3) {
      stop("Expected a 3d array of 'prob_strain'")
    }
    if (dim(prob_strain)[[2]] != length(pars)) {
      stop(sprintf(
        "Expected 2nd dim of 'prob_strain' to have length %d, following 'pars'",
        length(pars)))
    }
    if (dim(prob_strain)[[3]] != length(step)) {
      stop(sprintf(
        "Expected 3rd dim of 'prob_strain' to have length %d, following 'step'",
        length(step)))
    }
  }

  calculate_rt_one_trajectory <- function(i) {
    if (initial_step_from_parameters) {
      step[[1L]] <- pars[[i]]$initial_step
    }
    if (is.null(prob_strain)) {
      rt_1 <- calculate_Rt(step, S[, i, ], pars[[i]])
    } else {
      rt_1 <- calculate_Rt(step, S[, i, ], pars[[i]], prob_strain[, i, ])
    }
    rt_1
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

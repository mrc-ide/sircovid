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
##'   `prob_strain[1, j]` and `prob_strain[n_groups + 1, j]` should give, for
##'   the j^th time step, the probabilities that new infections
##'   in group 1 are of strains 1 and 2 respectively.
##'   The default is `NULL`, but it must
##'   be specified if there is more than one strain
##'
##' @param type A character vector of possible Rt types to
##'   compute. Can be any or all of `eff_Rt_all`, `eff_Rt_general`,
##'   `Rt_all` and `Rt_general`
##'
##' @param interpolate_every Spacing (in days) to use between interpolated
##'   points
##'
##' @param interpolate_critical_dates Optional vector of critical sircovid
##'   dates to use when interpolating. Interpolation will be done in
##'   blocks between the first time step of these dates, with each block
##'   starting on the first time step of the date given.
##'   So if you give a `interpolate_critical_dates` of `c(20,
##'   50)` then blocks *start* on first time step of days 20 and 50,
##'   i.e.: `[1, 20)`, `[20, 50)`, `[50, end]`.
##'
##' @param interpolate_min The minimum number of steps to include
##'   within a block. If there are fewer points than this then all
##'   points are used (i.e., no interpolation is done) or
##'   `interpolate_every` is reduced until at least this many points
##'   were used. This can be used to specify a lower bound on the
##'   error of small regions. If Rt is small it won't matter that
##'   much. You do need to specify something though or interpolation
##'   will not happen, and do not use less than 3 as we use spline
##'   interpolation and that will not work with fewer than 3 points.
##'
##' @param eigen_method The eigenvalue method to use (passed to
##'   [eigen1::eigen1] as `method`)
##'
##' @return A list with elements `step`, `beta`, and any of the `type`
##'   values specified above.
##'
##' @export
carehomes_Rt <- function(step, S, p, prob_strain = NULL,
                         type = NULL, interpolate_every = NULL,
                         interpolate_critical_dates = NULL,
                         interpolate_min = NULL,
                         eigen_method = "power_iteration") {
  all_types <- c("eff_Rt_all", "eff_Rt_general", "Rt_all", "Rt_general")
  if (is.null(type)) {
    type <- all_types
  } else {
    err <- setdiff(type, all_types)
    if (length(err) > 0) {
      stop(sprintf("Unknown R type %s, must match %s",
                   paste(squote(err), collapse = ", "),
                   paste(squote(all_types), collapse = ", ")))
    }
  }

  if (!is.null(interpolate_every)) {
    interpolate_critical_index <- match(interpolate_critical_dates / p$dt, step)
    step_index_split <- interpolate_grid_critical_x(seq_along(step),
                                                    interpolate_every,
                                                    interpolate_critical_index,
                                                    interpolate_min)
    step_index <- unlist(step_index_split)
    ret <- carehomes_Rt(step[step_index], S[, step_index, drop = FALSE], p,
                        prob_strain, type)
    if (!is.null(interpolate_every)) {
      ret[type] <- lapply(ret[type], interpolate_grid_expand_y,
                          step_index_split)
    }
    ## Also need to update these
    ret$step <- step
    ret$date <- step * p$dt
    ret$beta <- sircovid_parameters_beta_expand(step, p$beta_step)
    return(ret)
  }

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

  n_vacc_classes <- ncol(p$rel_susceptibility)
  n_strains <- length(p$strain_transmission)

  ### here mean_duration accounts for relative infectivity of
  ### different infection / vaccination stages
  beta <- sircovid_parameters_beta_expand(step, p$beta_step)

  ages <- seq_len(p$n_age_groups)
  ch <- seq(p$n_age_groups + 1L, p$n_groups)
  if (n_vacc_classes > 1L) {
    ch <- c(outer(ch, (seq_len(n_vacc_classes) - 1L) * p$n_groups, "+"))
  }

  ## We only need to do this section for each different beta value,
  ## and then it's still a pretty straightfoward scaling; we've
  ## applied beta to all *but* a fraction of the matrix (i.e., in the matrix
  ##
  ##   A B
  ##   C D
  ##
  ## Everything but D is scaled by beta
  ##
  ## When we have more than one vaccination group we replicate this to
  ## get block structure
  ##
  ##   A B A B
  ##   C D C D
  ##   A B A B
  ##   C D C D
  ##
  ## And scale all of A, B, C (not D) by beta
  m <- block_expand(unname(p$m), n_vacc_classes)
  mt <- m %o% beta
  mt[ch, ch, ] <- m[ch, ch]

  n_time <- length(step)
  prob_strain_mat <- array(prob_strain, c(p$n_groups, n_strains, n_time))
  if (any(is.na(prob_strain))) {
    ret <- list(step = step,
                date = step * p$dt,
                beta = beta)
    ret[type] <- list(rep(NA_real_, length(step)))
    return(ret)
  }

  mean_duration <-
    carehomes_Rt_mean_duration_weighted_by_infectivity(step, p,
                                                       prob_strain_mat)

  compute_ngm <- function(S) {
    len <- p$n_groups * n_vacc_classes
    Sw <- S * c(p$rel_susceptibility)
    mt * vapply(seq_len(n_time), function(t)
      tcrossprod(c(mean_duration[, , t]), Sw[, t]),
      matrix(0, len, len))
  }

  ## NOTE the signs on the exponents here is different! This gives
  ## good performance and reasonable accuracy to the point where
  ## this calculation is small in the profile.
  eigen <- function(m) {
    eigen1::eigen1(m, max_iterations = 1e5, tolerance = 1e-6,
                   method = eigen_method)
  }

  ret <- list(step = step,
              date = step * p$dt,
              beta = beta)

  if (any(c("eff_Rt_all", "eff_Rt_general") %in% type)) {
    ngm <- compute_ngm(S)
    if ("eff_Rt_all" %in% type) {
      ret$eff_Rt_all <- eigen(ngm)
    }
    if ("eff_Rt_general" %in% type) {
      ret$eff_Rt_general <- eigen(ngm[-ch, -ch, ])
    }
  }

  if (any(c("Rt_all", "Rt_general") %in% type)) {
    N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
    N_tot_all_vacc_groups <- N_tot_non_vacc
    if (n_vacc_classes > 1) {
      for (i in 2:n_vacc_classes) {
        N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                       0 * N_tot_non_vacc)
      }
    }
    ngm <- compute_ngm(N_tot_all_vacc_groups)
    if ("Rt_all" %in% type) {
      ret$Rt_all <- eigen(ngm)
    }
    if ("Rt_general" %in% type) {
      ret$Rt_general <- eigen(ngm[-ch, -ch, ])
    }
  }

  ret
}


## Here we expect 'S' in order:
##
##   state x sample x step
##
## We expect 'pars' to be a list along sample (or a shared parameter set)
## We expect 'step' to be a vector along step

##' Compute "Rt" for a set of simulated trajectories (e.g., the result
##' of the `$iterate()` method of [carehomes], [mcstate::pmcmc()] or
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
##' @inheritParams carehomes_Rt
##'
##' @return As for [carehomes_Rt()], except that every element is a
##'   matrix, not a vector.
##'
##' @export
carehomes_Rt_trajectories <- function(step, S, pars, prob_strain = NULL,
                                      initial_step_from_parameters = TRUE,
                                      shared_parameters = NULL,
                                      type = NULL,
                                      interpolate_every = NULL,
                                      interpolate_critical_dates = NULL,
                                      interpolate_min = NULL,
                                      eigen_method = "power_iteration") {
  calculate_Rt_trajectories(
    calculate_Rt = carehomes_Rt, step = step,
    S = S, pars = pars,
    prob_strain = prob_strain,
    initial_step_from_parameters = initial_step_from_parameters,
    shared_parameters = shared_parameters,
    type = type,
    interpolate_every = interpolate_every,
    interpolate_critical_dates = interpolate_critical_dates,
    interpolate_min = interpolate_min,
    eigen_method = eigen_method)
}


carehomes_Rt_mean_duration_weighted_by_infectivity <- function(step, pars,
                                                               strain_mat) {
  dt <- pars$dt

  matricise <- function(vect, n_col) {
    matrix(rep(vect, n_col), ncol = n_col, byrow = FALSE)
  }

  n_vacc_classes <- ncol(pars$rel_susceptibility)

  n_groups <- pars$n_groups

  n_time_steps <-
    length(sircovid_parameters_beta_expand(step, pars$p_H_step))

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

  ## just a helper function to calculate denominator after converting gammas
  ## to array
  dims <- c(dim(p_C), length(pars$gamma_A))
  n_groups <- dims[1L]
  n_vax <- dims[2L]
  n_time_steps <- dims[3L]
  n_strains <- dims[4L]

  denom <- function(x) array(rep(1 - exp(-dt * x), each = prod(dims[1:3])),
                             dims)

  mean_duration_I_A <- array(pars$I_A_transmission * (1 - p_C) *
    pars$k_A, dims) / denom(pars$gamma_A)

  mean_duration_I_P <- array(pars$I_P_transmission * p_C *
    pars$k_P, dims) / denom(pars$gamma_P)

  mean_duration_I_C_1 <- array(pars$I_C_1_transmission * p_C *
    pars$k_C_1, dims) / denom(pars$gamma_C_1)

  mean_duration_I_C_2 <- array(pars$I_C_2_transmission * p_C *
    pars$k_C_2, dims) / denom(pars$gamma_C_2)

  mean_duration_G_D <- array(pars$G_D_transmission * p_C * p_H *
    p_G_D * pars$k_G_D / (1 - exp(-dt * pars$gamma_G_D)), dims)

  mean_duration_hosp <- array(pars$hosp_transmission * (
    prob_H_R * pars$k_H_R / (1 - exp(-dt * pars$gamma_H_R)) +
      prob_H_D * pars$k_H_D / (1 - exp(-dt * pars$gamma_H_D)) +
      (prob_ICU_W_R + prob_ICU_W_D + prob_ICU_D) * pars$k_ICU_pre /
        (1 - exp(-dt * pars$gamma_ICU_pre))), dims)

  mean_duration_icu <- array(pars$ICU_transmission * (
    prob_ICU_W_R * pars$k_ICU_W_R / (1 - exp(-dt * pars$gamma_ICU_W_R)) +
      prob_ICU_W_D * pars$k_ICU_W_D / (1 - exp(-dt * pars$gamma_ICU_W_D)) +
      prob_ICU_D * pars$k_ICU_D / (1 - exp(-dt * pars$gamma_ICU_D))), dims)
  ## safety check
  stopifnot(identical(mean_duration_icu[, , , 1],
                      mean_duration_icu[, , , n_strains]))

  mean_duration <- mean_duration_I_A + mean_duration_I_P +
    mean_duration_I_C_1 + mean_duration_I_C_2 + mean_duration_G_D +
    mean_duration_hosp + mean_duration_icu

  ## Account for different infectivity levels depending on vaccination stage

  mean_duration <- mean_duration *
    outer(outer(pars$rel_infectivity, rep(1, n_time_steps)), rep(1, n_strains))

  ## Multiply by dt to convert from time steps to days
  mean_duration <- dt * mean_duration

  if (n_strains > 1L) {
    weighted_strain_multiplier <- vapply(seq_len(n_time_steps), function(t)
      matrix(strain_mat[,, t] %*% pars$strain_transmission,
             pars$n_groups, n_vax),
      matrix(0, pars$n_groups, n_vax))
    mean_duration <- mean_duration * outer(weighted_strain_multiplier,
                                           rep(1, n_strains))

    mean_duration <- vapply(seq_len(n_time_steps), function(t)
      vapply(seq_len(n_vax), function(j)
        rowSums(strain_mat[, , t] * mean_duration[, j, t, ]),
      numeric(dims[1L])),
      matrix(0, pars$n_groups, n_vax))
  } else {
    dim(mean_duration) <- dim(mean_duration)[1:3]
  }

  mean_duration
}


## This part holds over all possible Rt calculations, so I've factored
## it out here; when we implement this for the basic model this will
## remain unchanged.  However, I am leaving it in this
## carehomes-specific file until we do add a new model or port it.
calculate_Rt_trajectories <- function(calculate_Rt, step, S, pars, prob_strain,
                                      initial_step_from_parameters,
                                      shared_parameters, type, ...) {
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
      rt_1 <- calculate_Rt(step, S[, i, ], pars[[i]], type = type, ...)
    } else {
      rt_1 <- calculate_Rt(step, S[, i, ], pars[[i]], prob_strain[, i, ],
                           type = type, ...)
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

carehomes_Rt2 <- function(step, S, p, prob_strain = NULL,
                          type = NULL, interpolate_every = NULL,
                          interpolate_critical_dates = NULL,
                          interpolate_min = NULL) {
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
  mean_duration <- carehomes_Rt_mean_duration_weighted_by_infectivity2(step, p)

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

  ## TODO: do this rehape above?
  n_time <- length(step)
  prob_strain_mat <- array(prob_strain, c(p$n_groups, n_strains, n_time))
  if (any(is.na(prob_strain))) {
    stop("NA value in prob_strain - implement this")
  }

  ## TODO: if this is a timesink we can certainly do this with some
  ## bookkeeping trick especially as this works out to be a weighted
  ## sum.
  if (n_strains > 1L) {
    weighted_strain_multiplier <- vapply(seq_len(n_time), function(t)
      matrix(prob_strain_mat[, , t] %*% p$strain_transmission,
             p$n_groups, n_vacc_classes),
      matrix(0, p$n_groups, n_vacc_classes))
    mean_duration <- mean_duration * weighted_strain_multiplier
  }

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
    eigen1::eigen1(m, max_iterations = 1e5, tolerance = 1e-4)
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


carehomes_Rt_mean_duration_weighted_by_infectivity2 <- function(step, pars) {
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
  vacc_prog_before_infectious <- function(vacc_prog_rate) {
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

  ## mean_duration[i, j, k] represents mean duration at step k of age
  ## group i leaving the E compartment in vaccine stage j, we need to
  ## output for leaving the S compartment in vaccine stage j, so we
  ## calculate this here. If vaccine_progression_rate_base does not
  ## vary between age groups we lump them together here which saves a
  ## matrix inversion.
  if (n_vacc_classes > 1) {
    pr <- matrix_index(pars$vaccine_progression_rate_base)
    V <- vapply(pr$unique, function(i)
      vacc_prog_before_infectious(pr$value[i, ]),
      array(0, c(n_vacc_classes, n_vacc_classes)))[, , pr$index]

    out <- array(0, dim(mean_duration))
    for (i in seq_len(pars$n_groups)) {
      out[i, , ] <- V[, , i] %*% mean_duration[i, , ]
    }
  } else {
    out <- mean_duration
  }

  out
}

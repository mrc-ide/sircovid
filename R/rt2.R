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
  mean_duration <- carehomes_Rt_mean_duration_weighted_by_infectivity(step, p)

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
  n_time <- length(steps)
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

  N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
  N_tot_all_vacc_groups <- N_tot_non_vacc
  if (n_vacc_classes > 1) {
    for (i in 2:n_vacc_classes) {
      N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                     0 * N_tot_non_vacc)
    }
  }

  len <- p$n_groups * n_vacc_classes

  f <- function(S, drop_carehomes) {
    ## We weight S by relative susceptibility (only has an effect with
    ## vaccination); this is the only S used below. Because of the
    ## 'c()', these align nicely.
    Sw <- S * c(p$rel_susceptibility)

    ngm <- mt * vapply(seq_len(n_time), function(t)
      tcrossprod(c(mean_duration[, , t]), Sw[, t]),
      matrix(0, len, len))

    if (drop_carehomes) {
      ngm <- ngm[-ch, -ch, ]
    }

    ## NOTE the signs on the exponents here is different! This gives
    ## good performance and reasonable accuracy to the point where
    ## this calculation appears to vanish from the profile!
    eigen1::eigen1(ngm, max_iterations = 1e5, tolerance = 1e-4)
  }

  opts <- list(eff_Rt_all = list(pop = S, general = FALSE),
               eff_Rt_general = list(pop = S, general = TRUE),
               Rt_all = list(pop = N_tot_all_vacc_groups, general = FALSE),
               Rt_general = list(pop = N_tot_all_vacc_groups, general = TRUE))

  ret <- list(step = step,
              date = step * p$dt,
              beta = beta)

  ret[type] <- lapply(opts[type], function(x) f(x$pop, x$general))

  ret
}

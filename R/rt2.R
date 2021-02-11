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
  if (n_vacc_classes != 1) {
    stop("vaccination not supported")
  }
  if (length(p$strain_transmission) != 1) {
    stop("variants not supported")
  }

  ### here mean_duration accounts for relative infectivity of
  ### different infection / vaccination stages
  beta <- sircovid_parameters_beta_expand(step, p$beta_step)
  mean_duration <- carehomes_Rt_mean_duration_weighted_by_infectivity(step, p)

  ages <- seq_len(p$n_age_groups)
  ch <- seq(p$n_age_groups + 1L, p$n_groups) # 18:19

  ## We only need to do this section for each different beta value,
  ## and then it's still a pretty straightfoward scaling; we've
  ## applied beta to all *but* a fraction of the matrix (i.e., in the matrix
  ##
  ##   A B
  ##   C D
  ##
  ## Everything but D is scaled by beta; in the case of an expanded
  ## matrix this requires a bit more work
  mt <- unname(p$m) %o% beta
  mt[ch, ch, ] <- p$m[ch, ch]

  N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
  N_tot_all_vacc_groups <- N_tot_non_vacc
  if (n_vacc_classes > 1) {
    for (i in 2:n_vacc_classes) {
      N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                     0 * N_tot_non_vacc)
    }
  }

  f <- function(S, drop_carehomes) {
    ## We weight S by relative susceptibility (only has an effect with
    ## variants); this is the only S used below.
    S_weighted <- S * c(p$rel_susceptibility)

    ## TODO: in the presence of strains there is some work to do to get
    ## mean duration correct.
    ngm <- mt * vapply(seq_along(steps), function(t)
      tcrossprod(mean_duration[, , t], S_weighted[, t]),
      matrix(0, 19, 19))

    if (drop_carehomes) {
      ## Lots of ways of generating this
      i <- rep(rep(c(TRUE, FALSE), c(length(ages), length(ch))),
               length.out = nrow(ngm))
      ngm <- ngm[i, i, ]
    }

    ## NOTE the signs on the exponents here is different! This gives
    ## good performance and reasonable accuracy to the point where
    ## this calculation appears to vanish from the profile!
    eigen1::eigen1(ngm, max_iterations = 1e5, tolerance = 1e-5)
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

##' Compute "Rt" for a single simulated trajectory and parameter set.
##'
##' @title Compute "Rt"
##'
##' @param time A vector of time steps that the model was run over
##'
##' @param S A (n groups x n vaccine classes) x steps matrix of "S"
##'   compartment counts
##'
##' @param p A [lancelot_parameters()] object
##'
##' @param prob_strain A (n groups x n strains) x n time steps matrix of
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
##' @param R A (n groups x n strains x n vaccine classes) x time steps matrix of
##'   "R" compartment counts, required for multi-strain models.
##'
##' @param weight_Rt If `TRUE` then computes the weighted average
##'  of the Rt for all strains, otherwise all calculations are returned with
##'  an additional dimension to index each strain.
##'
##' @param keep_strains_Rt Additional argument for when `weight_Rt` is `TRUE`
##'  (has no impact otherwise). If `TRUE`, then the Rt for each strain is
##'  returned along with the weighted average, otherwise just the weighted
##'  average is returned. When `TRUE`, the dimension indexing these lists
##'  strains first, and then the weighted average.
##'
##' @return A list with elements `time`, `beta`, and any of the `type`
##'   values specified above.
##'
##' @export
lancelot_Rt <- function(time, S, p, prob_strain = NULL,
                        type = NULL, interpolate_every = NULL,
                        interpolate_critical_dates = NULL,
                        interpolate_min = NULL,
                        eigen_method = "power_iteration", R = NULL,
                        weight_Rt = FALSE, keep_strains_Rt = FALSE) {

  if (sum(p$hosp_transmission, p$ICU_transmission, p$G_D_transmission) > 0) {
    stop("Cannot currently compute Rt if any of 'hosp_transmission',
    'ICU_transmission' or 'G_D_transmission' are non-zero")
  }

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

  n_strains <- length(p$strain_transmission)
  ## move prob_strain check up here and make NULL if not needed to shortcut
  ##  checks and calculations
  if (n_strains > 1) {
    if (is.null(prob_strain) || is.null(R)) {
      stop("Expected prob_strain and R input because there is more than one
            strain")
      ## deal with the prob_strain NA case at the beginning before any variables
      ##  are modified
    } else if (!is.null(prob_strain) && any(is.na(prob_strain)) && weight_Rt) {
      which_nna <- !vlapply(seq(ncol(prob_strain)),
                            function(i) any(is.na(prob_strain[, i])))

      ## calculate Rt by first ignoring NA
      if (any(which_nna)) {
        ret <- lancelot_Rt(
          time[which_nna], S[, which_nna, drop = FALSE], p,
          prob_strain[, which_nna, drop = FALSE], type, interpolate_every,
          interpolate_critical_dates, interpolate_min,
          eigen_method, R[, which_nna, drop = FALSE], weight_Rt,
          keep_strains_Rt)
      } else {
        ret <- vector("list", 3 + length(type))
        names(ret) <- c("time", "date", "beta", type)
      }

      ## replace reduced time,date,beta with full values (no NA here)
      ret$time <- time
      ret$date <- time * p$dt
      ret$beta <- beta <- sircovid_parameters_expand_step(time, p$beta_step)

      ## restore full length Rt with NAs when prob_strain is NA
      for (i in grep("Rt_", names(ret))) {
        base <- rep(NA, length(time))
        base[which_nna] <- ret[[i]]
        ret[[i]] <- base
      }

      return(ret)
    }
  } else {
    prob_strain <- array(1, length(time))
    R <- NULL
  }

  if (!is.null(interpolate_every)) {
    interpolate_critical_index <- match(interpolate_critical_dates / p$dt,
                                        time)
    # remove NA values and ensure vector order is decreasing
    interpolate_critical_index <-
      interpolate_critical_index[!is.na(interpolate_critical_index)]
    interpolate_critical_index <-
      interpolate_critical_index[order(interpolate_critical_index,
                                       decreasing = FALSE)]

    time_index_split <- interpolate_grid_critical_x(seq_along(time),
                                                    interpolate_every,
                                                    interpolate_critical_index,
                                                    interpolate_min)
    time_index <- unlist(time_index_split)
    ret <- lancelot_Rt(time[time_index], S[, time_index, drop = FALSE], p,
                       prob_strain[, time_index, drop = FALSE], type,
                       R = R[, time_index, drop = FALSE],
                       weight_Rt = weight_Rt, keep_strains_Rt = keep_strains_Rt)
    if (!is.null(interpolate_every)) {
      ret[type] <- lapply(ret[type], interpolate_grid_expand_y,
                          time_index_split)
    }
    ## Also need to update these
    ret$time <- time
    ret$date <- time * p$dt
    ret$beta <- sircovid_parameters_expand_step(time, p$beta_step)
    return(ret)
  }

  if (nrow(S) != nlayer(p$rel_susceptibility) * nrow(p$m)) {
    stop(sprintf(
      "Expected 'S' to have %d rows = %d groups x %d vaccine classes",
      p$n_groups * nlayer(p$rel_susceptibility),
      p$n_groups,
      nlayer(p$rel_susceptibility)))
  }
  if (ncol(S) != length(time)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'time'",
                 length(time)))
  }
  if (!is.null(R)) {
    if (nrow(R) != nlayer(p$rel_susceptibility) * nrow(p$m) * p$n_strains_R) {
      stop(sprintf(
        "Expected 'R' to have %d rows = %d groups x %d strains_R x %d vaccine
          classes",
        p$n_groups * nlayer(p$rel_susceptibility) * p$n_strains_R,
        p$n_groups, p$n_strains_R, nlayer(p$rel_susceptibility)))
    }
    if (ncol(R) != length(time)) {
      stop(sprintf("Expected 'R' to have %d columns, following 'time'",
                   length(time)))
    }
  }

  if (n_strains == 1) {
    n_real_strains <- 1
  } else {
    ## unmirror pseudo-strains value (with safety checks)
    p <- unmirror_pars(p)
    n_real_strains <- 2

    if (!is.matrix(prob_strain)) {
      stop(sprintf(
        "Expected a %d strains x %d time steps matrix for 'prob_strain'",
        n_real_strains, length(time)))
    }

    if (nrow(prob_strain) != n_real_strains) {
      stop(sprintf(
        "Expected 'prob_strain' to have %d rows, following number of strains",
        n_real_strains))
    }
    if (ncol(prob_strain) != length(time)) {
      stop(sprintf(
        "Expected 'prob_strain' to have %d columns, following 'time'",
        length(time)))
    }
  }

  n_vacc_classes <- nlayer(p$rel_susceptibility)

  ### here mean_duration accounts for relative infectivity of
  ### different infection / vaccination stages
  beta <- sircovid_parameters_expand_step(time, p$beta_step)

  ages <- seq_len(p$n_age_groups)
  if (p$has_carehomes == 1) {
    ch <- seq(p$n_age_groups + 1L, p$n_groups)
    if (n_vacc_classes > 1L) {
      ch <- c(outer(ch, (seq_len(n_vacc_classes) - 1L) * p$n_groups, "+"))
    }
  }

  ## We only need to do this section for each different beta value,
  ## and then it's still a pretty straightforward scaling; we've
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
  if (p$has_carehomes == 1) {
    mt[ch, ch, ] <- m[ch, ch]
  }

  n_time <- length(time)

  mean_duration <-
    lancelot_Rt_mean_duration_weighted_by_infectivity(time, p)

  compute_ngm <- function(x, S, rel_sus, R = 0, rel_sus_strain = 1) {
    n_groups <- p$n_groups
    len <- n_groups * n_vacc_classes
    Sw <- (S + R * rel_sus_strain) * c(rel_sus)
    mt * vapply(seq_len(n_time), function(t)
      tcrossprod(c(x[, , t]), Sw[, t]),
      matrix(0, len, len))
  }

  ## NOTE the signs on the exponents here is different! This gives
  ## good performance and reasonable accuracy to the point where
  ## this calculation is small in the profile.
  eigen <- function(m) {
    eigen1::eigen1(m, max_iterations = 1e5, tolerance = 1e-6,
                   method = eigen_method)
  }

  ret <- list(time = time,
              date = time * p$dt,
              beta = beta)

  n_groups <- nrow(p$m)

  ngm_computer <- function(x, effective) {
    vapply(seq(n_real_strains), function(i) {
      md <- mean_duration[, i, , , drop = FALSE]
      dim(md) <- dim(md)[c(1, 3, 4)]

      if (n_real_strains == 1) {
        compute_ngm(md, x, p$rel_susceptibility)
      } else {
        N <- n_groups * n_vacc_classes
        ## in the new model set-up it can only be 5, can make a variable
        ##  if this changes in the future
        RR <- array(R, c(n_groups, p$n_strains_R, n_vacc_classes, ncol(R)))

        if (i == 1) {
          if (effective) {
            ## get R5 (as they're susceptible to E1)
            RR <- matrix(RR[, 5, , ], nrow = n_groups * n_vacc_classes)
          } else {
            RR <- 0
          }
        } else {
          if (effective) {
            ## get R1, R4 and R5 (as they're susceptible to E2)
            RR <- matrix(RR[, 1, , ], nrow = n_groups * n_vacc_classes) +
              matrix(RR[, 4, , ], nrow = n_groups * n_vacc_classes) +
              matrix(RR[, 5, , ], nrow = n_groups * n_vacc_classes)
          } else {
            RR <- 0
          }
        }
        ## We are calculating the NGM for strain i, so we need the cross
        ## immunity of strain 3 - i (2 if i = 1, 1 if i = 2) against strain i
        compute_ngm(md, x, p$rel_susceptibility[, i, ], RR,
                    1 - p$cross_immunity[3 - i])
      }
    }, array(0, c(n_groups * n_vacc_classes, n_groups * n_vacc_classes,
                  n_time)))
  }

  if (any(c("eff_Rt_all", "eff_Rt_general") %in% type)) {
    ngm <- ngm_computer(S, effective = TRUE)

    if ("eff_Rt_all" %in% type) {
      ret$eff_Rt_all <-
        vapply(seq(n_real_strains), function(i) eigen(ngm[, , , i]),
               numeric(dim(ngm)[[3L]]))
    }
    if ("eff_Rt_general" %in% type) {
      if (p$has_carehomes == 1) {
        ret$eff_Rt_general <-
          vapply(seq(n_real_strains), function(i) eigen(ngm[-ch, -ch, , i]),
                 numeric(dim(ngm)[[3L]]))
      } else {
        ret$eff_Rt_general <-
          vapply(seq(n_real_strains), function(i) eigen(ngm[, , , i]),
                 numeric(dim(ngm)[[3L]]))
      }

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

    ngm <- ngm_computer(N_tot_all_vacc_groups, effective = FALSE)

    if ("Rt_all" %in% type) {
      ret$Rt_all <-
        vapply(seq(n_real_strains), function(i) eigen(ngm[, , , i]),
               numeric(dim(ngm)[[3L]]))
    }
    if ("Rt_general" %in% type) {
      if (p$has_carehomes == 1) {
        ret$Rt_general <-
          vapply(seq(n_real_strains), function(i) eigen(ngm[-ch, -ch, , i]),
                 numeric(dim(ngm)[[3L]]))
      } else {
        ret$Rt_general <-
          vapply(seq(n_real_strains), function(i) eigen(ngm[, , , i]),
                 numeric(dim(ngm)[[3L]]))
      }
    }
  }

  ## ensure backwards compatibility by dropping columns for single_strain and
  ## separating classes
  class(ret) <- c("Rt")
  is_single <- (inherits(last(ret), c("matrix", "array")) &&
                  is.null(ncol(last(ret)))) ||
    (!inherits(last(ret), c("matrix", "array")) &&
       length(last(ret)) == 1)


  if (is_single) {
    ## adding 'nocov' as this is a safety check that should never be hit
    class(ret) <- c("single_strain", class(ret)) # nocov
  } else if (isTRUE(ncol(last(ret)) == 1)) {
    ret[intersect(all_types, names(ret))] <-
      lapply(ret[intersect(all_types, names(ret))], drop)
    class(ret) <- c("single_strain", class(ret))
  } else {
    if (weight_Rt) {
      ret <- wtmean_Rt(ret, prob_strain, keep_strains_Rt)
      if (keep_strains_Rt) {
        class(ret) <- c("multi_strain_weighted", class(ret))
      } else {
        ## treat multi strain as single once weighted
        class(ret) <- c("single_strain", class(ret))
      }

    } else {
      class(ret) <- c("multi_strain", class(ret))
    }
  }

  ret
}


## Here we expect 'S' in order:
##
##   state x sample x time
##
## We expect 'pars' to be a list along sample (or a shared parameter set)
## We expect 'time' to be a vector along time

##' Compute "Rt" for a set of simulated trajectories (e.g., the result
##' of the `$iterate()` method of [lancelot], [mcstate::pmcmc()] or
##' [mcstate::pmcmc_predict()]. The trajectories may or may not share
##' parameters.
##'
##' @title Compute Rt for a set of trajectories
##'
##' @param time A vector of time steps
##'
##' @param S A 3d ((n groups x n vaccine classes) x n trajectories x n time
##'   steps) array of "S" compartment counts
##'
##' @param pars Either a single [lancelot_parameters()] object
##'   (shared parameters) or an unnamed list of
##'   [lancelot_parameters()] objects, the same length as `ncol(S)`.
##'
##' @param prob_strain A 3d ((n groups x n strains) x n trajectories x n time
##'   steps) array of "prob_strain" model outputs. Default is `NULL`, but it
##'   must be specified if there is more than one strain.
##'
##' @param initial_time_from_parameters If `TRUE`, then `time[[1]]` is
##'   replaced by the value of `initial_time` from the parameters.
##'   This is usually what you want. (From sircovid 0.12.13 this
##'   parameter means "initial time is zero" and will probably be
##'   updated in a future version).
##'
##' @param shared_parameters Should `pars` be treated as a single
##'   shared list? Leave as `NULL` to detect automatically, set to
##'   `TRUE` or `FALSE` to force it to be interpreted one way or the
##'   other which may give more easily interpretable error messages.
##'
##' @param R A 3d ((n groups x n strains x n vaccine classes) x
##'   n trajectories x n time steps) array of "R" compartment counts, required
##'   for multi-strain models.
##'
##' @inheritParams lancelot_Rt
##'
##' @return As for [lancelot_Rt()], except that every element is a
##'   matrix, not a vector.
##'
##' @export
lancelot_Rt_trajectories <- function(time, S, pars, prob_strain = NULL,
                                     initial_time_from_parameters = TRUE,
                                     shared_parameters = NULL,
                                     type = NULL,
                                     interpolate_every = NULL,
                                     interpolate_critical_dates = NULL,
                                     interpolate_min = NULL,
                                     eigen_method = "power_iteration",
                                     R = NULL, weight_Rt = FALSE,
                                     keep_strains_Rt = FALSE) {
  calculate_Rt_trajectories(
    calculate_Rt = lancelot_Rt, time = time,
    S = S, pars = pars,
    prob_strain = prob_strain,
    initial_time_from_parameters = initial_time_from_parameters,
    shared_parameters = shared_parameters,
    type = type,
    interpolate_every = interpolate_every,
    interpolate_critical_dates = interpolate_critical_dates,
    interpolate_min = interpolate_min,
    eigen_method = eigen_method,
    R = R,
    weight_Rt = weight_Rt,
    keep_strains_Rt = keep_strains_Rt)
}


lancelot_Rt_mean_duration_weighted_by_infectivity <- function(time, pars) {

  dt <- pars$dt
  n_time_steps <-
    length(sircovid_parameters_expand_step(time, pars$p_H_step))

  p_C <- combine_steps_groups(
    time, pars$n_groups, n_time_steps,
    n_strains = length(pars$strain_transmission),
    n_vacc_classes = pars$n_vacc_classes, p_step = pars$p_C_step,
    rel_p = pars$rel_p_sympt,
    strain_rel_p = pars$strain_rel_p_sympt
  )

  ## Compute mean duration (in time steps) of each stage of infection,
  ## weighed by probability of going through that stage
  ## and by relative infectivity of that stage

  ## Note the mean duration (in time steps) of a compartment for
  ## a discretised Erlang(k, gamma) is k / (1 - exp(dt * gamma))
  calculate_mean <- function(transmission, prob, name) {
    gamma_step <-
      sircovid_parameters_expand_step(time,
                                      pars[[paste0("gamma_", name, "_step")]])
    rel_gamma <- pars[[paste0("rel_gamma_", name)]]
    k <- pars[[paste0("k_", name)]]
    gamma <- aperm(outer(outer(gamma_step, rel_gamma),
                         array(1, c(pars$n_groups, pars$n_vacc_classes))),
                   c(3, 2, 4, 1))
    transmission * k * prob / stats::pexp(gamma, dt)
  }

  mean_duration_I_A <- calculate_mean(pars$I_A_transmission, (1 - p_C), "A")
  mean_duration_I_P <- calculate_mean(pars$I_P_transmission, p_C, "P")
  mean_duration_I_C_1 <- calculate_mean(pars$I_C_1_transmission, p_C, "C_1")
  mean_duration_I_C_2 <- calculate_mean(pars$I_C_2_transmission, p_C, "C_2")

  mean_duration <- mean_duration_I_A + mean_duration_I_P +
    mean_duration_I_C_1 + mean_duration_I_C_2

  ## Account for different infectivity levels depending on vaccination stage
  mean_duration <- mean_duration * outer(pars$rel_infectivity,
                                         rep(1, n_time_steps))

  ## Account for different infectivity levels depending on strain
  ##  reorder to match dimensions of mean_duration and strain_transmission then
  ##  revert back
  mean_duration <-
    aperm(aperm(mean_duration, c(2, 1, 3, 4)) * pars$strain_transmission,
          c(2, 1, 3, 4))

  ## Multiply by dt to convert from time steps to days
  mean_duration <- dt * mean_duration

  mean_duration
}

## This part holds over all possible Rt calculations, so I've factored
## it out here; when we implement this for the basic model this will
## remain unchanged.  However, I am leaving it in this
## lancelot-specific file until we do add a new model or port it.
calculate_Rt_trajectories <- function(calculate_Rt, time, S, pars, prob_strain,
                                      initial_time_from_parameters,
                                      shared_parameters, type, R = NULL, ...) {
  if (length(dim(S)) != 3) {
    stop("Expected a 3d array of 'S'")
  }

  if (!is.null(R) && length(dim(R)) != 3) {
    stop("Expected a 3d array of 'R'")
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

  if (dim(S)[[3]] != length(time)) {
    stop(sprintf(
      "Expected 3rd dimension of 'S' to have length %d, following 'time'",
      length(time)))
  }

  if (!is.null(R) && any(dim(R)[2:3] != dim(S)[2:3])) {
    stop(sprintf(
      "Expected 2nd and 3rd dimension of 'R' to be the same as 'S' (%d x %d)'",
      ncol(S), dim(S)[[3]]))
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
    if (dim(prob_strain)[[3]] != length(time)) {
      stop(sprintf(
        "Expected 3rd dim of 'prob_strain' to have length %d, following 'time'",
        length(time)))
    }
  }

  calculate_rt_one_trajectory <- function(i) {
    if (initial_time_from_parameters) {
      ## TODO: Ed (or someone else) this has been probably not ideal
      ## since we moved to seeding as this is *always* zero now!
      ## Similar problem in the ifr calculation.
      ##
      ## The current formulation will be backward compatible and leave
      ## tests passing until we fix this properly.
      time[[1L]] <- pars[[i]]$initial_time %||% 0
    }
    if (is.null(prob_strain)) {
      rt_1 <- calculate_Rt(time, S[, i, ], pars[[i]], type = type,
                           R = R[, i, ], ...)
    } else {
      rt_1 <- calculate_Rt(time, S[, i, ], pars[[i]], prob_strain[, i, ],
                           type = type, R = R[, i, ], ...)
    }
    rt_1
  }

  res <- lapply(seq_along(pars), calculate_rt_one_trajectory)

  ## These are stored in a list-of-lists and we convert to a
  ## list-of-matrices here
  collect <- function(nm) {
    if (nm %in% c("time", "date", "beta")) {
      matrix(unlist(lapply(res, "[[", nm)), length(time), length(res))
    } else {
      array(unlist(lapply(res, "[[", nm)),
            c(length(time), ncol(res[[1]][[nm]]), length(res)))
    }
  }
  nms <- names(res[[1]])
  ret <- set_names(lapply(nms, collect), nms)

  ## ensure backwards compatibility by dropping columns for single_strain and
  ## separating classes
  all_types <- c("eff_Rt_all", "eff_Rt_general", "Rt_all", "Rt_general")
  if (length(dim(ret[[length(ret)]])) < 3) {
    ret[intersect(all_types, names(ret))] <-
      lapply(ret[intersect(all_types, names(ret))], drop)
    class(ret) <- c("single_strain", "Rt_trajectories", "Rt")
  } else if (dim(ret[[length(ret)]])[2] < 3) {
    class(ret) <- c("multi_strain", "Rt_trajectories", "Rt")
  } else {
    class(ret) <- c("multi_strain_weighted", "Rt_trajectories", "Rt")
  }

  ret
}


wtmean_Rt <- function(rt, prob_strain, keep_strains_Rt) {
  if (!inherits(rt, "Rt")) {
    stop("'rt' must inherit from class 'Rt")
  }

  rt_mean <- rt

  what <- names(rt)
  what <- what[!(what %in% c("time", "date", "beta"))]

  get_mean_rt <- function(r, prob_strain) {
    n_dim <- length(dim(prob_strain))
    strain_dim <- 2
    reshape_prob_strain <- aperm(prob_strain,
                                 c(n_dim, seq_len(n_dim)[-n_dim]))

    if (!inherits(r, c("matrix", "array"))) {
      r <- matrix(r, nrow = 1)
    }

    if (length(dim(r)) != length(dim(reshape_prob_strain)) ||
        !all(dim(r) == dim(reshape_prob_strain))) {
      stop(sprintf(
        "Expect elements of Rt to have dimensions: %d time steps x %d strains x
        %d particles", nrow(reshape_prob_strain), ncol(reshape_prob_strain),
        nlayer(reshape_prob_strain)))
    }
    ## catch the case when strain_transmission is 0 so r is NaN
    x <- r * reshape_prob_strain
    x[reshape_prob_strain == 0] <- 0
    res <- apply(x, seq_len(n_dim)[-strain_dim], sum)
    if (keep_strains_Rt) {
      res <- cbind(r, res)
    }
    res
  }

  rt_mean[what] <- lapply(what, function(i) get_mean_rt(rt[[i]], prob_strain))

  rt_mean

}

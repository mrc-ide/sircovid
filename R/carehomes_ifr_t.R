##' Compute "IFR_t" for a single simulated trajectory and parameter set.
##'
##' @title Compute "IFR_t"
##'
##' @param step A vector of steps that the model was run over
##'
##' @param S A (n groups x n vaccine classes) x steps matrix of
##'   "S" compartment counts
##'
##' @param I_weighted A (n groups x n strains x n vaccine classes) x steps
##'   matrix of "I_weighted" compartment counts
##'
##' @param p A [carehomes_parameters()] object
##'
##' @param type A character vector of possible Rt types to
##'   compute. Can be any or all of `IFR_t_all`, `IFR_t_general`,
##'   `IHR_t_all`, `IHR_t_general`, `IFR_t_all_no_vacc`,
##'   `IFR_t_general_no_vacc`, `IHR_t_all_no_vacc` and `IHR_t_general_no_vacc`,
##'
##' @param R A (n groups x n strains x n vaccine classes) x steps matrix of "R"
##'   compartment counts, required for multi-strain models.
##'
##'
##' @return A list with elements `step`, `date`, and any of the `type`
##'   values specified above.
##'
##' @importFrom stats weighted.mean
##'
##' @export
carehomes_ifr_t <- function(step, S, I_weighted, p, type = NULL, R = NULL) {
  all_types <- c("IFR_t_all", "IFR_t_general", "IHR_t_all", "IHR_t_general",
                 "IFR_t_all_no_vacc", "IFR_t_general_no_vacc",
                 "IHR_t_all_no_vacc", "IHR_t_general_no_vacc",
                 "ALOS", "ALOS_no_vacc")
  if (is.null(type)) {
    type <- all_types
  } else {
    err <- setdiff(type, all_types)
    if (length(err) > 0) {
      stop(sprintf("Unknown IFR/IHR type %s, must match %s",
                   paste(squote(err), collapse = ", "),
                   paste(squote(all_types), collapse = ", ")))
    }
  }

  n_groups <- p$n_groups
  n_strains <- p$n_strains
  n_vacc_classes <- p$n_vacc_classes

  if (n_strains > 1) {
    if (n_strains != 4) {
      stop("Multstrain IFR currently only works if n_strains is 4")
    }
    if (is.null(R)) {
      stop("Expected R input because there is more than one
            strain")
    }
    if (nrow(R) != n_groups * n_strains * n_vacc_classes) {
      stop(sprintf(
        "Expected 'R' to have %d rows = %d groups x %d strains x %d
          vaccine classes",
        n_groups * n_strains * n_vacc_classes,
        n_groups, n_strains, n_vacc_classes))
    }
    if (ncol(R) != length(step)) {
      stop(sprintf("Expected 'R' to have %d columns, following 'step'",
                   length(step)))
    }
  }

  if (nrow(S) != n_groups * n_vacc_classes) {
    stop(sprintf(
      "Expected 'S' to have %d rows = %d groups x %d vacc classes",
      n_groups * n_vacc_classes, n_groups, n_vacc_classes))
  }
  if (ncol(S) != length(step)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'step'",
                 length(step)))
  }
  if (nrow(I_weighted) != n_groups * n_strains * n_vacc_classes) {
    stop(sprintf(
      "Expected 'I_weighted' to have %d rows = %d groups x %d strains x %d
          vaccine classes",
      n_groups * n_strains * n_vacc_classes,
      n_groups, n_strains, n_vacc_classes))
  }
  if (ncol(I_weighted) != length(step)) {
    stop(sprintf("Expected 'I_weighted' to have %d columns, following 'step'",
                 length(step)))
  }

  if (n_strains == 1) {
    n_real_strains <- 1
  } else {
    ## unmirror pseudo-strains value (with safety checks)
    p <- unmirror_pars(p)
    n_real_strains <- 2
  }

  IFR_by_group_and_vacc_class <-
    carehomes_IFR_t_by_group_and_vacc_class(step, p)
  beta <- sircovid_parameters_expand_step(step, p$beta_step)

  calculate_expected_infections <- function(t, no_vacc) {
    m <- p$m
    ages <- seq_len(p$n_age_groups)
    ch <- seq(to = n_groups, length.out = 2)
    m[ages, ] <- beta[t] * m[ages, ]
    m[ch, ages] <- beta[t] * m[ch, ages]

    m_extended <- matrix(t(matrix(m, n_groups, n_groups * n_vacc_classes)),
                         n_groups * n_vacc_classes,
                         n_groups * n_vacc_classes,
                         byrow = TRUE)

    ## probability of infection in next time step by group and vaccine class
    rel_sus <- p$rel_susceptibility
    rel_inf <- p$rel_infectivity
    if (no_vacc) {
      rel_sus[, , ] <- 1
      rel_inf[, , ] <- 1
    }

    ## expected infections in next time step by group and vaccine class
    if (n_strains == 1) {
      ## force of infection by group and vaccine class (as a vector)
      foi <- c(rel_sus) * (m_extended %*% (I_weighted[, t] * c(rel_inf)))
      ## probability of infection by group and vaccine class
      prob_infection <- (1 - exp(-p$dt * foi))
      ## expected infections by group and vaccine class
      expected_infections <- S[, t] * prob_infection
    } else {
      ## need to sum over metastrains
      II <- array(I_weighted[, t], c(n_groups, n_strains, n_vacc_classes))
      II[, 1, ] <- II[, 1, ] + II[, 4, ]
      II[, 2, ] <- II[, 2, ] + II[, 3, ]
      II <- II[, -c(3, 4), , drop = FALSE]

      ## foi[i, j, k] is force of infection on a susceptible in group i/vaccine
      ## class k from strain j
      foi <- array(0, c(n_groups, n_real_strains, n_vacc_classes))
      for (i in 1:2) {
        foi_strain <- c(rel_sus[, i, ]) *
          (m_extended %*% (c(II[, i, ]) *
                             c(p$strain_transmission[i]) * c(rel_inf[, i, ])))
        foi[, i, ] <- array(foi_strain, c(n_groups, n_vacc_classes))
      }

      SS <- array(S[, t], c(n_groups, n_vacc_classes))
      RR <- array(R[, t], c(n_groups, n_strains, n_vacc_classes))
      ## total force of infection across strains
      total_foi <- apply(foi, c(1, 3), sum)
      ## calculate expected infections by group, strain and vaccine classes
      ## accounting for cross immunity
      e_inf <- array(0, c(n_groups, n_real_strains, n_vacc_classes))
      for (i in 1:2) {
        other_strain <- ifelse(i == 1, 2, 1)
        e_inf[, i, ] <- SS * (1 - exp(-p$dt * total_foi)) *
          foi[, i, ] / total_foi + RR[, other_strain, ] *
          (1 - exp(-p$dt * foi[, i, ] * (1 - p$cross_immunity[other_strain])))
      }
      expected_infections <- c(e_inf)
    }
    expected_infections
  }


  calculate_weighted_ratio <- function(t, expected_infections,
                                       drop_carehomes, no_vacc, type) {

    ## Care home workers (CHW) and residents (CHR) in last two rows
    ## and columns, remove for each vaccine class
    if (drop_carehomes) {
      i_CHR <- seq(n_groups, n_groups * n_real_strains * n_vacc_classes,
                   by = n_groups)
      i_CHW <- i_CHR - 1
      i_keep <-
        seq_len(n_groups * n_real_strains * n_vacc_classes)[-c(i_CHW, i_CHR)]
    } else {
      i_keep <- seq_len(n_groups * n_real_strains * n_vacc_classes)
    }

    if (no_vacc) {
      ## same IFR by group across all vaccine classes
      IFR_vec <- rep(c(IFR_by_group_and_vacc_class[[type]][, , 1, t]),
                     n_vacc_classes)
    } else {
      IFR_vec <- c(IFR_by_group_and_vacc_class[[type]][, , , t])
    }

    if (type != "ALOS") {
      out <- weighted.mean(IFR_vec[i_keep], expected_infections[i_keep, t])
    } else {
      if (no_vacc) {
        ## same IHR by group across all vaccine classes
        IHR_vec <- rep(c(IFR_by_group_and_vacc_class[["IHR"]][, , 1, t]),
                       n_vacc_classes)
      } else {
        IHR_vec <- c(IFR_by_group_and_vacc_class[["IHR"]][, , , t])
      }
      out <- weighted.mean(IFR_vec[i_keep],
                           expected_infections[i_keep, t] * IHR_vec[i_keep])
    }
    out
  }


  t <- seq_along(step)
  expected_infections_vacc <- vapply(t, calculate_expected_infections,
                                     numeric(n_groups * n_real_strains *
                                               n_vacc_classes),
                                     no_vacc = FALSE)
  expected_infections_no_vacc <- vapply(t, calculate_expected_infections,
                                        numeric(n_groups * n_real_strains *
                                                  n_vacc_classes),
                                        no_vacc = TRUE)

  opts <- list(IFR_t_all = list(e_inf = expected_infections_vacc,
                                general = FALSE, no_vacc = FALSE,
                                type = "IFR"),
               IFR_t_general = list(e_inf = expected_infections_vacc,
                                    general = TRUE, no_vacc = FALSE,
                                    type = "IFR"),
               IHR_t_all = list(e_inf = expected_infections_vacc,
                                general = FALSE, no_vacc = FALSE,
                                type = "IHR"),
               IHR_t_general = list(e_inf = expected_infections_vacc,
                                    general = TRUE, no_vacc = FALSE,
                                    type = "IHR"),
               IFR_t_all_no_vacc = list(e_inf = expected_infections_no_vacc,
                                        general = FALSE, no_vacc = TRUE,
                                        type = "IFR"),
               IFR_t_general_no_vacc = list(e_inf = expected_infections_no_vacc,
                                            general = TRUE, no_vacc = TRUE,
                                            type = "IFR"),
               IHR_t_all_no_vacc = list(e_inf = expected_infections_no_vacc,
                                        general = FALSE, no_vacc = TRUE,
                                        type = "IHR"),
               IHR_t_general_no_vacc = list(e_inf = expected_infections_no_vacc,
                                            general = TRUE, no_vacc = TRUE,
                                            type = "IHR"),
               ALOS = list(e_inf = expected_infections_vacc,
                           general = FALSE, no_vacc = FALSE, type = "ALOS"),
               ALOS_no_vacc = list(e_inf = expected_infections_no_vacc,
                                   general = FALSE, no_vacc = TRUE,
                                   type = "ALOS"))

  ret <- list(step = step,
              date = step * p$dt)
  ret[type] <- lapply(opts[type], function(x)
    vnapply(t, calculate_weighted_ratio, x$e_inf,
            drop_carehomes = x$general,
            no_vacc = x$no_vacc, type = x$type))
  ret

}

## Here we expect 'S' and 'I_weighted' in order:
##
##   state x sample x step
##
## We expect 'pars' to be a list along sample (or a shared parameter set)
## We expect 'step' to be a vector along step

##' Compute "IFR_t" for a set of simulated trajectories (e.g., the result
##' of the `$iterate()` method of [carehomes], [mcstate::pmcmc()] or
##' [mcstate::pmcmc_predict()]. The trajectories may or may not share
##' parameters.
##'
##' @title Compute IFR_t for a set of trajectories
##'
##' @param step A vector of steps
##'
##' @param S A 3d ((n groups x n vaccine classes) x n trajectories
##'   x n steps) array of "S" compartment counts
##'
##' @param I_weighted A 3d ((n groups x n strains x n vaccine classes) x
##'   n trajectories x n steps) array of "I_weighted" compartment counts
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
##' @param R A 3d ((n groups x n strains x n vaccine classes) x
##'   n trajectories x n steps) array of "R" compartment counts, required for
##'   multi-strain models.
##'
##' @inheritParams carehomes_ifr_t
##'
##' @return As for [carehomes_ifr_t()], except that every element is a
##'   matrix, not a vector.
##'
##' @export
carehomes_ifr_t_trajectories <- function(step, S, I_weighted, pars,
                                      initial_step_from_parameters = TRUE,
                                      shared_parameters = NULL, type = NULL,
                                      R = NULL) {
  calculate_ifr_t_trajectories(carehomes_ifr_t, step, S, I_weighted, pars,
                            initial_step_from_parameters, shared_parameters,
                            type, R)
}


carehomes_IFR_t_by_group_and_vacc_class <- function(step, pars) {

  probs <- compute_pathway_probabilities(
    step = step,
    pars = pars,
    n_time_steps = length(sircovid_parameters_expand_step(step,
                                                          pars$p_H_step)),
    n_strains = length(pars$strain_transmission),
    n_vacc_classes = nlayer(pars$rel_susceptibility))

  p_C <- probs$p_C
  p_H <- probs$p_H
  p_ICU <- probs$p_ICU
  p_ICU_D <- probs$p_ICU_D
  p_H_D <- probs$p_H_D
  p_W_D <- probs$p_W_D
  p_G_D <- probs$p_G_D

  dt <- pars$dt

  ## Note the mean duration (in time steps) of a compartment for
  ## a discretised Erlang(k, gamma) is k / (1 - exp(dt * gamma))
  calculate_mean <- function(name) {
    gamma_step <-
      sircovid_parameters_expand_step(step,
                                      pars[[paste0("gamma_", name, "_step")]])
    rel_gamma <- pars[[paste0("rel_gamma_", name)]]
    k <- pars[[paste0("k_", name)]]
    gamma <- aperm(outer(outer(gamma_step, rel_gamma),
                         array(1, c(pars$n_groups, pars$n_vacc_classes))),
                   c(3, 2, 4, 1))
    k / stats::pexp(gamma, dt)
  }

  duration_ICU_pre <- calculate_mean("ICU_pre")
  duration_ICU_D <- calculate_mean("ICU_D")
  duration_ICU_W_D <- calculate_mean("ICU_W_D")
  duration_ICU_W_R <- calculate_mean("ICU_W_R")
  duration_H_D <- calculate_mean("H_D")
  duration_H_R <- calculate_mean("H_R")
  duration_W_D <- calculate_mean("W_D")
  duration_W_R <- calculate_mean("W_R")

  IHR <- p_C * p_H * (1 - p_G_D) * 100
  IFR <- p_C * p_H * p_G_D * 100 +
    IHR * (p_ICU * (p_ICU_D + (1 - p_ICU_D) * p_W_D) +
             (1 - p_ICU) * p_H_D)

  ALOS <- (1 - p_ICU) * (p_H_D * duration_H_D + (1 - p_H_D) * duration_H_R) +
    p_ICU * (duration_ICU_pre + (1 - p_ICU_D) * duration_ICU_D +
               p_W_D * (duration_ICU_W_D + duration_W_D) +
               (1 - p_W_D) * (duration_ICU_W_R + duration_W_R))
  ALOS <- dt * ALOS

  out <- list(IFR = IFR,
              IHR = IHR,
              ALOS = ALOS)

  out
}


## This part holds over all possible IFR_t calculations, so I've factored
## it out here; when we implement this for the basic model this will
## remain unchanged.  However, I am leaving it in this
## carehomes-specific file until we do add a new model or port it.
calculate_ifr_t_trajectories <- function(calculate_ifr_t, step, S, I_weighted,
                                         pars, initial_step_from_parameters,
                                         shared_parameters, type, R = NULL) {
  if (length(dim(S)) != 3) {
    stop("Expected a 3d array of 'S'")
  }
  if (length(dim(I_weighted)) != 3) {
    stop("Expected a 3d array of 'I_weighted'")
  }
  if (!is.null(R) && length(dim(R)) != 3) {
    stop("Expected a 3d array of 'R'")
  }

  shared_parameters <- shared_parameters %||% !is.null(names(pars))
  if (shared_parameters) {
    if (is.null(names(pars))) {
      stop("If using shared parameters, expected a named list for 'pars'")
    }
    if (dim(I_weighted)[[2]] != dim(S)[[2]]) {
      stop("Expected 'S' and 'I_weighted' to have same length of 2nd dim")
    }
    if (!is.null(R) && dim(R)[[2]] != dim(S)[[2]]) {
      stop("Expected 'S' and 'R' to have same length of 2nd dim")
    }
    pars <- rep(list(pars), ncol(S))
  } else {
    if (!is.null(names(pars))) {
      stop("If not using shared parameters, expected a unnamed list for 'pars'")
    }
    if (length(pars) != ncol(S)) {
      stop(sprintf(
        "Expected 2nd dim of 'S' to have length %d, given 'pars'",
        length(pars)))
    }
    if (length(pars) != ncol(I_weighted)) {
      stop(sprintf(
        "Expected 2nd dim of 'I_weighted' to have length %d, given 'pars'",
        length(pars)))
    }
    if (!is.null(R) && length(pars) != ncol(R)) {
      stop(sprintf(
        "Expected 2nd dim of 'R' to have length %d, given 'pars'",
        length(pars)))
    }
  }

  if (dim(S)[[3]] != length(step)) {
    stop(sprintf(
      "Expected 3rd dim of 'S' to have length %d, given 'step'",
      length(step)))
  }
  if (dim(I_weighted)[[3]] != length(step)) {
    stop(sprintf(
      "Expected 3rd dim of 'I_weighted' to have length %d, given 'step'",
      length(step)))
  }
  if (!is.null(R) && dim(R)[[3]] != length(step)) {
    stop(sprintf(
      "Expected 3rd dim of 'R' to have length %d, given 'step'",
      length(step)))
  }

  calculate_ifr_t_one_trajectory <- function(i) {
    if (initial_step_from_parameters) {
      step[[1L]] <- pars[[i]]$initial_step
    }
    if (is.null(R)) {
      ifr_t_1 <- calculate_ifr_t(step, S[, i, ], I_weighted[, i, ], pars[[i]],
                                 type = type)
    } else {
      ifr_t_1 <- calculate_ifr_t(step, S[, i, ], I_weighted[, i, ], pars[[i]],
                                 type = type, R = R[, i, ])
    }

    ifr_t_1
  }

  res <- lapply(seq_along(pars), calculate_ifr_t_one_trajectory)

  ## These are stored in a list-of-lists and we convert to a
  ## list-of-matrices here
  collect <- function(nm) {
    matrix(unlist(lapply(res, "[[", nm)), length(step), length(res))
  }
  nms <- names(res[[1]])
  ret <- set_names(lapply(nms, collect), nms)
  class(ret) <- "IFR_t_trajectories"
  ret
}

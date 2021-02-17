##' Compute "IFR_t" for a single simulated trajectory and parameter set.
##'
##' @title Compute "IFR_t"
##'
##' @param step A vector of steps that the model was run over
##'
##' @param S A (n groups x n vaccine classes) x steps matrix of
##'   "S" compartment counts
##'
##' @param I_weighted A (n groups x n vaccine classes) x steps matrix of
##'   "I_weighted" compartment counts
##'
##' @param p A [carehomes_parameters()] object
##'
##' @param type A character vector of possible Rt types to
##'   compute. Can be any or all of `IFR_t_all`, `IFR_t_general`,
##'   `IHR_t_all`, `IHR_t_general`, `IFR_t_all_no_vacc`,
##'   `IFR_t_general_no_vacc`, `IHR_t_all_no_vacc` and `IHR_t_general_no_vacc`,
##'
##' @return A list with elements `step`, `date`, and any of the `type`
##'   values specified above.
##'
##' @importFrom stats weighted.mean
##'
##' @export
carehomes_ifr_t <- function(step, S, I_weighted, p, type = NULL) {

  all_types <- c("IFR_t_all", "IFR_t_general", "IHR_t_all", "IHR_t_general",
                 "IFR_t_all_no_vacc", "IFR_t_general_no_vacc",
                 "IHR_t_all_no_vacc", "IHR_t_general_no_vacc")
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


  if (nrow(S) != ncol(p$rel_susceptibility) * p$n_groups) {
    stop(sprintf(
      "Expected 'S' to have %d rows = %d groups x %d vacc classes",
      p$n_groups * ncol(p$rel_susceptibility),
      p$n_groups,
      ncol(p$rel_susceptibility)))
  }
  if (ncol(S) != length(step)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'step'",
                 length(step)))
  }
  if (nrow(I_weighted) != ncol(p$rel_susceptibility) * p$n_groups) {
    stop(sprintf(
      "Expected 'I_weighted' to have %d rows = %d groups x %d vacc classes",
      p$n_groups * ncol(p$rel_susceptibility),
      p$n_groups,
      ncol(p$rel_susceptibility)))
  }
  if (ncol(I_weighted) != length(step)) {
    stop(sprintf("Expected 'I_weighted' to have %d columns, following 'step'",
                 length(step)))
  }

  IFR_by_group_and_vacc_class <-
    carehomes_IFR_t_by_group_and_vacc_class(step, p)
  beta <- sircovid_parameters_beta_expand(step, p$beta_step)
  n_vacc_classes <- ncol(p$rel_susceptibility)


  calculate_expected_infections <- function(t, no_vacc) {
    m <- p$m
    ages <- seq_len(p$n_age_groups)
    ch <- seq(to = p$n_groups, length.out = 2)
    m[ages, ] <- beta[t] * m[ages, ]
    m[ch, ages] <- beta[t] * m[ch, ages]

    m_extended <- matrix(t(matrix(m, p$n_groups, p$n_groups * n_vacc_classes)),
                         p$n_groups * n_vacc_classes,
                         p$n_groups * n_vacc_classes,
                         byrow = TRUE)

    ## probability of infection in next time step by group and vaccine class
    if (no_vacc) {
      ## exclude vaccine effects
      prob_infection <- (1 - exp(-p$dt * (m_extended %*% I_weighted[, t])))
    } else {
      ## include vaccine effects
      prob_infection <- (1 - exp(-p$dt * c(p$rel_susceptibility) *
                     (m_extended %*% (I_weighted[, t] * c(p$rel_infectivity)))))
    }

    ## expected infections in next time step by group and vaccine class
    expected_infections <-  S[, t] * prob_infection

    expected_infections
  }


  calculate_weighted_ratio <- function(t, expected_infections,
                                       drop_carehomes, no_vacc, type) {

    ## Care home workers (CHW) and residents (CHR) in last two rows
    ## and columns, remove for each vaccine class
    if (drop_carehomes) {
      i_CHR <- seq(p$n_groups, dim(S)[1], by = p$n_groups)
      i_CHW <- i_CHR - 1
      i_keep <- seq_len(dim(S)[1])[-c(i_CHW, i_CHR)]
    } else {
      i_keep <- seq_len(dim(S)[1])
    }

    if (no_vacc) {
      ## same IFR by group across all vaccine classes
      IFR_vec <- rep(c(IFR_by_group_and_vacc_class[[type]][, 1, t]),
                     n_vacc_classes)
    } else {
      IFR_vec <- c(IFR_by_group_and_vacc_class[[type]][, , t])
    }

    weighted.mean(IFR_vec[i_keep], expected_infections[i_keep, t])
  }


  t <- seq_along(step)
  expected_infections_vacc <- vapply(t, calculate_expected_infections,
                                     numeric(dim(S)[1]), no_vacc = FALSE)
  expected_infections_no_vacc <- vapply(t, calculate_expected_infections,
                                        numeric(dim(S)[1]), no_vacc = TRUE)

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
                                            type = "IHR"))

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
##' @param I_weighted A 3d ((n groups x n vaccine classes) x n trajectories
##'   x n steps) array of "I_weighted" compartment counts
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
##' @inheritParams carehomes_ifr_t
##'
##' @return As for [carehomes_ifr_t()], except that every element is a
##'   matrix, not a vector.
##'
##' @export
carehomes_ifr_t_trajectories <- function(step, S, I_weighted, pars,
                                      initial_step_from_parameters = TRUE,
                                      shared_parameters = NULL, type = NULL) {
  calculate_ifr_t_trajectories(carehomes_ifr_t, step, S, I_weighted, pars,
                            initial_step_from_parameters, shared_parameters,
                            type)
}


carehomes_IFR_t_by_group_and_vacc_class <- function(step, pars) {
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

  IHR <- p_C * p_H * (1 - p_G_D) * 100
  IFR <- p_C * p_H * p_G_D * 100 +
    IHR * (p_ICU * (p_ICU_D + (1 - p_ICU_D) * p_W_D) +
             (1 - p_ICU) * p_H_D)

  out <- list(IFR = IFR,
              IHR = IHR)
  
  out
}


## This part holds over all possible IFR_t calculations, so I've factored
## it out here; when we implement this for the basic model this will
## remain unchanged.  However, I am leaving it in this
## carehomes-specific file until we do add a new model or port it.
calculate_ifr_t_trajectories <- function(calculate_ifr_t, step, S, I_weighted,
                                         pars, initial_step_from_parameters,
                                         shared_parameters, type) {
  if (length(dim(S)) != 3) {
    stop("Expected a 3d array of 'S'")
  }
  if (length(dim(I_weighted)) != 3) {
    stop("Expected a 3d array of 'I_weighted'")
  }

  shared_parameters <- shared_parameters %||% !is.null(names(pars))
  if (shared_parameters) {
    if (is.null(names(pars))) {
      stop("If using shared parameters, expected a named list for 'pars'")
    }
    if (dim(I_weighted)[[2]] != dim(S)[[2]]) {
      stop("Expected 'S' and 'I_weighted' to have same length of 2nd dim")
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

  calculate_ifr_t_one_trajectory <- function(i) {
    if (initial_step_from_parameters) {
      step[[1L]] <- pars[[i]]$initial_step
    }
    ifr_t_1 <- calculate_ifr_t(step, S[, i, ], I_weighted[, i, ], pars[[i]],
                               type = type)
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

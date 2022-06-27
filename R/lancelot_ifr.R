##' Compute "IFR" excluding immunity.
##'
##' @title Compute "IFR" excluding immunity
##'
##' @param step A vector of steps to calculate IFR at
##'
##' @param pars An unnamed list of [lancelot_parameters()] objects
##'
##' @export
lancelot_ifr_excl_immunity <- function(step, pars) {

  if (length(dim(step)) > 1) {
    stop("Expected 'step' to be a vector")
  }
  if (length(unique(unlist(lapply(pars, function(x) x$n_strains)))) > 1) {
    stop("All parameter sets must have the same number of strains")
  }
  if (length(unique(unlist(lapply(pars, function(x) x$n_age_groups)))) > 1) {
    stop("All parameter sets must have the same number of age groups")
  }

  if (pars[[1]]$n_strains == 1) {
    n_real_strains <- 1
  } else {
    n_real_strains <- 2
  }
  n_age_groups <- pars[[1]]$n_age_groups
  n_pars <- length(pars)

  gen_pop <- seq_len(n_age_groups)

  compute_ifr_weighting_strain <- function(p, i_strain) {
    if (sum(p$hosp_transmission, p$ICU_transmission, p$G_D_transmission) > 0) {
      stop("Cannot currently compute IFR if any of 'hosp_transmission',
           'ICU_transmission' or 'G_D_transmission' are non-zero")
    }

    mean_duration <-
      lancelot_Rt_mean_duration_weighted_by_infectivity(step, p)

    ifr_weight_step <- function(t) {
      ## We want unvaccinated only so set 1 in 3rd dimension
      ngm <- p$m[gen_pop, gen_pop] *
        tcrossprod(mean_duration[gen_pop, i_strain, 1, t], p$N_tot[gen_pop])

      weighting <- rep(0, p$n_age_groups)
      ## We need to exclude empty groups from the eigenvector calculation
      ## They will be weighted 0
      k <- which(p$N_tot[gen_pop] != 0)
      weighting[k] <- eigen(ngm[k, k])$vectors[, 1]

      assert_real(weighting)
    }

    vapply(seq_along(step), ifr_weight_step, numeric(n_age_groups))
  }

  compute_ifr_weighting <- function(p) {
    out <- vapply(seq_len(n_real_strains),
                  function(i) compute_ifr_weighting_strain(p, i),
                  array(0, c(n_age_groups, length(step))))
    out <- aperm(out, c(1, 3, 2))
    out
  }

  ifr_weighting <-
    vapply(pars, compute_ifr_weighting,
           array(0, c(n_age_groups, n_real_strains, length(step))))

  ifr_unweighted <-
    lapply(pars, function(p) lancelot_ifr_by_group_strain_vacc_class(step, p))

  weight_ifr <- function(i_strain, i_step, i_pars, type) {
    if (type == "HFR") {
      ## for HFR, need to weight further by IHR
      w <- ifr_weighting[, i_strain, i_step, i_pars] *
        ifr_unweighted[[i_pars]]$IHR[gen_pop, i_strain, 1, i_step]
    } else {
      w <- ifr_weighting[, i_strain, i_step, i_pars]
    }
    weighted.mean(
      ifr_unweighted[[i_pars]][[type]][gen_pop, i_strain, 1, i_step], w)
  }

  calc_type <- function(type) {
    ## This will output an object of dim: n_steps x n_real_strains x n_pars
    out <-
      vapply(seq_len(n_pars),
             function(i_pars)
               vapply(seq_len(n_real_strains),
                      function(i_strain)
                        vnapply(seq_along(step),
                                function(i_step)
                                  weight_ifr(i_strain, i_step, i_pars, type)),
                    numeric(length(step))),
            array(0, c(length(step), n_real_strains)))

    if (sum(c(length(step), n_real_strains, n_pars) == 1) >= 2) {
      out <- array(out, c(length(step), n_real_strains, n_pars))
    }
    out

  }

  ret <- list(IFR = calc_type("IFR"),
              IHR = calc_type("IHR"),
              HFR = calc_type("HFR"),
              step = step)

  ret
}

lancelot_ifr_by_group_strain_vacc_class <- function(step, pars) {

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

  IHR <- p_C * p_H * (1 - p_G_D)
  HFR <- p_ICU * (p_ICU_D + (1 - p_ICU_D) * p_W_D) +
    (1 - p_ICU) * p_H_D
  IFR <- p_C * p_H * p_G_D + IHR * HFR

  out <- list(IFR = IFR,
              IHR = IHR,
              HFR = HFR)


  out
}


compute_pathway_probabilities <- function(step, pars, n_time_steps, n_strains,
                                          n_vacc_classes) {

  i <- seq_len(n_strains)

  out <- list()
  out$p_C <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_C_step, pars$rel_p_sympt[, i, , drop = FALSE],
    pars$strain_rel_p_sympt)
  out$p_H <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_H_step, pars$rel_p_hosp_if_sympt[, i, , drop = FALSE],
    pars$strain_rel_p_hosp_if_sympt)
  out$p_ICU <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_ICU_step, pars$rel_p_ICU[, i, , drop = FALSE],
    pars$strain_rel_p_icu)
  out$p_ICU_D <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_ICU_D_step, pars$rel_p_ICU_D[, i, , drop = FALSE],
    pars$strain_rel_p_ICU_D)
  out$p_H_D <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_H_D_step, pars$rel_p_H_D[, i, , drop = FALSE],
    pars$strain_rel_p_H_D)
  out$p_W_D <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_W_D_step, pars$rel_p_W_D[, i, , drop = FALSE],
    pars$strain_rel_p_W_D)
  out$p_G_D <- combine_steps_groups(
    step, pars$n_groups, n_time_steps, n_strains, n_vacc_classes,
    pars$p_G_D_step, pars$rel_p_G_D[, i, , drop = FALSE],
    pars$strain_rel_p_G_D)

  out
}

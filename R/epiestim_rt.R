##' @title Performs draws from the generation time distribution
##'
##' @param p A single [carehomes_parameters()] object.
##'
##' @param n The number of draws to perform
##'
##' @return a vector of `n` values of the generation time
##'
##' @importFrom distcrete distcrete
##' @importFrom stats runif rpois
gt_sample <- function(p, n = 1000) {

  ## Note: the above does not account for the impact of
  ## vaccination on the generation time.

  if (length(unique(p$p_C)) > 1)  {
    stop("gt_sample does not allow p_C to vary by age")
  }

  if (p$I_C_2_transmission > 0) {
    stop("gt_sample does not allow transmission from I_C_2")
  }
  if (p$k_A > 1 || p$k_P  > 1 || p$k_C_1 > 1) {
    stop("gt_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1")
  }
  if (p$I_P_transmission != 1 || p$I_C_1_transmission  != 1) {
    stop("gt_sample does not allow
    I_P_transmission !=1 or I_C_1_transmission != 1")
  }

  draw_from <- function(n, k, gamma) {
    discrete_gamma <-
      distcrete::distcrete("gamma", p$dt, shape = k, rate = gamma, w = 1)
    discrete_gamma$r(n)
  }

  sampled_gt <- numeric(0)
  i <- 0

  while (i < n) {
    ## Draw exposed period
    sample_E <- draw_from(1, p$k_E, p$gamma_E)
    ## Draw symptomatic vs asymptomatic path
    path_C <- runif(1) <= p$p_C[[1]]
    if (!path_C) {
      sample_I <- draw_from(1, p$k_A, p$gamma_A)
      infectivity <- p$I_A_transmission
    } else {
      sample_I <- draw_from(1, p$k_P, p$gamma_P) +
        draw_from(1, p$k_C_1, p$gamma_C_1)
      infectivity <- 1
    }
    n_secondary_cases <- rpois(1, lambda = infectivity * sample_I)
    tmp_gt <- seq(sample_E + p$dt, sample_E + sample_I, p$dt)
    sampled_gt <- c(sampled_gt,
                    sample(tmp_gt, n_secondary_cases, replace = TRUE))
    i <- length(sampled_gt)
  }

  sampled_gt

}

##' @importFrom graphics hist
gt_distr <- function(p, n = 1000) {
  sampled_gt <- gt_sample(p, n)
  ## the discretisation below allows having a zero on the first day
  ## which is required by EpiEstim
  ## this is because the GT is computed from the sum of E and I
  ## and each of those is at least 0.25 days
  ## hence the GT is at least 0.5 day
  hist(sampled_gt, breaks = seq(0, ceiling(max(sampled_gt) + 1), 1) - 0.51,
       plot = FALSE)$density
}


## We expect 'pars' to be a list along sample
## We expect 'step' to be a vector along step

##' Compute "Rt" using EpiEtim for a set of simulated trajectories (e.g., the
##' result of the `$iterate()` method of [carehomes], [mcstate::pmcmc()] or
##' [mcstate::pmcmc_predict()]. The trajectories should share
##' parameters.
##'
##' @title Compute Rt using EpiEstim for a set of trajectories
##'
##' @param step A vector of steps
##'
##' @param incidence A matrix (n trajectories x n steps) of incidence counts
##'
##' @param p A single [carehomes_parameters()] object.
##'
##' @param sliding_window_ndays An integer giving the length of the sliding
##' window on which Rt will be estimated
##'
##' @param mean_prior The mean prior for Rt
##'
##' @param sd_prior The standard deviation of the prior for Rt
##'
##' @param n_GT An integer giving the number of generation times to be drawn to
##' construct the discrete distribution of the generation time
##'
##' @param n_R An integer giving the number of Rt values to sample from for each
##' incidence trajectory. These will then be aggregated across all incidence
##' trajectories.
##'
##' @param save_all A boolean determining whether to save all samples of Rt
##' estimated or only a summary
##'
##' @return A list with elements
##' `t_start` (vector of first days of the sliding windows over which Rt is
##' estimated),
##' `t_end` (vector of last days of the sliding windows over which Rt is
##' estimated),
##' `Rt` a matrix (only present if `save_all = TRUE`) containing for each
##' sliding window (each row in the matrix) a sample of n_R * nrow(inc) values
##' of Rt for that sliding window (columns of the matrix)
##' `Rt_summary` a matrix containing for each sliding window (each row in the
##' matrix) the 2.5%, 50%, 97.5% quantiles and the mean of Rt for that sliding
##' window (columns of the matrix)
##'
##' @export
##' @importFrom EpiEstim make_config estimate_R
##' @importFrom stats quantile rgamma
carehomes_EpiEstim_Rt_trajectories <- function(step, incidence, p,
                                               sliding_window_ndays = 7,
                                               mean_prior = 1,
                                               sd_prior = 1,
                                               n_GT = 10000,
                                               n_R = 1000,
                                               save_all = TRUE) {
  gt_distr <- gt_distr(p = p, n = n_GT)

  T <- ncol(incidence)
  np <- nrow(incidence)
  t_start <- seq(2, T - sliding_window_ndays + 1)
  t_end <- seq(sliding_window_ndays + 1, T)

  # generate the config based on incidence for first parameter set
  config <- EpiEstim::make_config(incid = incidence[1, ],
                                  t_start = t_start,
                                  t_end = t_end,
                                  method = "non_parametric_SI",
                                  si_distr = gt_distr,
                                  mean_prior = mean_prior,
                                  std_prior = sd_prior)

  R_sample <- matrix(NA, length(t_start), n_R * np)

  for (i in seq_len(np)) {
    ## This function gives warning when precision in estimates is not good
    ## in case users may over-interpret central estimate without considering
    ## uncertainty.
    ## Warning is not needed in our context as we will fully capture uncertainty
    R_i <- suppressWarnings(EpiEstim::estimate_R(incid = incidence[i, ],
                                                 config = config)$R)

    R_i_shape_scale <- gamma_mucv2shapescale(mu = R_i[["Mean(R)"]],
                                     cv = R_i[["Std(R)"]] / R_i[["Mean(R)"]])
    f_sample_R <- function(e) {
      if (!is.na(R_i_shape_scale$shape[e])) {
        ret <- rgamma(n_R, shape = R_i_shape_scale$shape[e],
                      scale = R_i_shape_scale$scale[e])
      } else {
        ret <- rep(NA, n_R)
      }
    }

    R_sample[, n_R * (i - 1) + seq_len(n_R)] <-
      t(vapply(seq_len(length(t_start)), f_sample_R, numeric(n_R)))
  }

  summary_R <- apply(R_sample, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
  mean_R <- apply(R_sample, 1, mean, na.rm = TRUE)
  summary_R <- t(rbind(summary_R, mean_R))

  time_start <- step[t_start] * p$dt
  time_end <- step[t_end] * p$dt

  if (save_all) {
    ret <- list(t_start = time_start,
                t_end = time_end,
                Rt = R_sample,
                Rt_summary = summary_R)
  } else {
    ret <- list(t_start = time_start,
                t_end = time_end,
                Rt_summary = summary_R)
  }

}
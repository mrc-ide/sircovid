test_that("Can calculate EpiEstim Rt (I)", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)

  np <- 20L
  mod <- carehomes$new(p, 0, np, seed = 2L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_cum_inf <- mod$info()$index$cum_infections
  index_S <- mod$info()$index$S
  index <- c(index_cum_inf, index_S)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  inc <- apply(y[1, , ], 1, diff)
  inc <- t(rbind(rep(0, np), inc))
  S <- y[-1, , ]

  rt_EpiEstim <- carehomes_rt_trajectories_epiestim(steps, inc, p,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  #### General patterns
  ## dimension of Rt should be 3 rows and length(steps) - 1 cols
  ## the -1 because EpiEstim only start estimation at 2nd time step
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(4, length(steps) - 1))
  ## Rt at the end should be < 1
  expect_vector_lt(rt_EpiEstim$Rt_summary[, ncol(rt_EpiEstim$Rt_summary)], 1)
  ## because of the priors with mean 1 we would expect EpiEstim Rt
  ## at the end to be higher than eff_Rt_all
  expect_true(
    last(rt$eff_Rt_all) <
      rt_EpiEstim$Rt_summary["mean_R", ncol(rt_EpiEstim$Rt_summary)])
  ## Rt at the start should be > 1
  ## (but there will be uncertainty so looking at the mean)
  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary["mean_R", ])))
  expect_true(rt_EpiEstim$Rt_summary["mean_R", first_non_NA_idx] > 1)
})

test_that("Can calculate EpiEstim Rt (II)", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)

  np <- 20L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_cum_inf <- mod$info()$index$cum_infections
  index_S <- mod$info()$index$S
  index <- c(index_cum_inf, index_S)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  inc <- apply(y[1, , ], 1, diff)
  inc <- t(rbind(rep(0, np), inc))
  S <- y[-1, , ]

  rt_EpiEstim <- carehomes_rt_trajectories_epiestim(steps, inc, p,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary["mean_R", ])))
  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", first_non_NA_idx]], 3.1,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", first_non_NA_idx]], 0.7,
               tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[["mean_R", ncol(rt_EpiEstim$Rt_summary)]], 0.2,
    tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[["97.5%", ncol(rt_EpiEstim$Rt_summary)]], 0.2,
    tolerance = .1)
})


test_that("Can calculate EpiEstim Rt with predefined GT (I)", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)
  np <- 20L
  mod <- carehomes$new(p, 0, np, seed = 2L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_cum_inf <- mod$info()$index$cum_infections
  index_S <- mod$info()$index$S
  index <- c(index_cum_inf, index_S)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  inc <- apply(y[1, , ], 1, diff)
  inc <- t(rbind(rep(0, np), inc))
  S <- y[-1, , ]

  gt_distr <- gt_distr(p, 10000)

  rt_EpiEstim <- carehomes_rt_trajectories_epiestim(steps, inc, p, gt_distr,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  #### General patterns
  ## dimension of Rt should be 3 rows and length(steps) - 1 cols
  ## the -1 because EpiEstim only start estimation at 2nd time step
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(4, length(steps) - 1))
  ## Rt at the end should be < 1
  expect_vector_lt(rt_EpiEstim$Rt_summary[, ncol(rt_EpiEstim$Rt_summary)], 1)
  ## because of the priors with mean 1 we would expect EpiEstim Rt
  ## at the end to be higher than eff_Rt_all
  expect_true(
    last(rt$eff_Rt_all) <
      rt_EpiEstim$Rt_summary["mean_R", ncol(rt_EpiEstim$Rt_summary)])
  ## Rt at the start should be > 1
  ## (but there will be uncertainty so looking at the mean)
  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary["mean_R", ])))
  expect_true(rt_EpiEstim$Rt_summary["mean_R", first_non_NA_idx] > 1)
})


test_that("Can calculate EpiEstim Rt with predefined GT (II)", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)
  np <- 20L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_cum_inf <- mod$info()$index$cum_infections
  index_S <- mod$info()$index$S
  index <- c(index_cum_inf, index_S)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  inc <- apply(y[1, , ], 1, diff)
  inc <- t(rbind(rep(0, np), inc))
  S <- y[-1, , ]

  gt_distr <- gt_distr(p, 10000)

  rt_EpiEstim <- carehomes_rt_trajectories_epiestim(steps, inc, p, gt_distr,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary["mean_R", ])))
  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", first_non_NA_idx]], 3.1,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", first_non_NA_idx]], 0.7,
               tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[["mean_R", ncol(rt_EpiEstim$Rt_summary)]], 0.2,
    tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[["97.5%", ncol(rt_EpiEstim$Rt_summary)]], 0.2,
    tolerance = .1)
})


test_that("carehomes_rt_trajectories_epiestim rejects invalid input", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)
  np <- 20L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_cum_inf <- mod$info()$index$cum_infections
  index_S <- mod$info()$index$S
  index <- c(index_cum_inf, index_S)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  inc <- apply(y[1, , ], 1, diff)
  inc <- t(rbind(rep(0, np), inc))
  S <- y[-1, , ]

  gt_distr <- gt_distr(p, 10000)
  gt_distr_mat <- matrix(rep(gt_distr, 2), nrow = 2, byrow = TRUE)
  gt_distr_not_zero <- gt_distr
  gt_distr_not_zero[1] <- 0.1

  gt_distr_not_sum_to_1 <- gt_distr * 2

  expect_error(carehomes_rt_trajectories_epiestim(steps, inc, p,
                                                    gt_distr_mat,
                                                    sliding_window_ndays = 1),
               "Expecting a vector for 'gt_distr'")

  expect_error(carehomes_rt_trajectories_epiestim(steps, inc, p,
                                                  gt_distr_not_zero,
                                                  sliding_window_ndays = 1),
               "The first value of 'gt_distr' should be zero.")

  expect_error(carehomes_rt_trajectories_epiestim(steps, inc, p,
                                                  gt_distr_not_sum_to_1,
                                                  sliding_window_ndays = 1),
               "'gt_distr' should sum to 1.")

})


test_that("Can calculate EpiEstim Rt when no transmission in carehomes", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)
  mean_N_tot <- mean(p$N_tot)
  mean_m <- mean(p$m)
  ## switch off transmission in carehomes
  p$m[, 18:19] <- 0
  p$m[18:19, ] <- 0
  ## make transmission homogeneous everywhere else
  mean(p$m)
  p$m[- (18:19), - (18:19)] <- mean_m
  ## same population size as well
  p$N_tot <- c(rep(mean_N_tot, 17), 0, 0)

  np <- 20L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_cum_inf <- mod$info()$index$cum_infections
  index_S <- mod$info()$index$S
  index <- c(index_cum_inf, index_S)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  inc <- apply(y[1, , ], 1, diff)
  inc <- t(rbind(rep(0, np), inc))
  S <- y[-1, , ]

  rt_EpiEstim <- carehomes_rt_trajectories_epiestim(steps, inc, p,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  #### General patterns
  ## dimension of Rt should be 3 rows and length(steps) - 1 cols
  ## the -1 because EpiEstim only start estimation at 2nd time step
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(4, length(steps) - 1))
  ## Rt at the end should be < 1
  expect_vector_lt(rt_EpiEstim$Rt_summary[, ncol(rt_EpiEstim$Rt_summary)], 1)
  ## because of the priors with mean 1 we would expect EpiEstim Rt
  ## at the end to be higher than eff_Rt_all
  expect_true(
    last(rt$eff_Rt_all) <
      rt_EpiEstim$Rt_summary["mean_R", ncol(rt_EpiEstim$Rt_summary)])
  ## Rt at the start should be > 1
  ## (but there will be uncertainty so looking at the mean)
  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary["mean_R", ])))
  expect_true(rt_EpiEstim$Rt_summary["mean_R", first_non_NA_idx] > 1)
  ## Rt at the start should be similar between methods:
  expect_true(
    rt$eff_Rt_all[1] > rt_EpiEstim$Rt_summary["2.5%", first_non_NA_idx])
  expect_true(
    rt$eff_Rt_all[1] < rt_EpiEstim$Rt_summary["97.5%", first_non_NA_idx])

  ## Rt stays constant for a while but precision improves for EpiEstim
  expect_true(rt$eff_Rt_all[1] > rt_EpiEstim$Rt_summary["2.5%", 20])
  expect_true(rt$eff_Rt_all[1] < rt_EpiEstim$Rt_summary["97.5%", 20])
  expect_true(abs(rt$eff_Rt_all[1] - rt_EpiEstim$Rt_summary["50%", 20]) < 0.2)

  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", 35]], 4.9,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", 35]], 4.8,
               tolerance = .1)
  expect_equal(rt$eff_Rt_all[35], 5.0,
               tolerance = .1)

  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", 45]], 4.2,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", 45]], 3.6,
               tolerance = .1)
  expect_equal(rt$eff_Rt_all[45], 4.57,
               tolerance = .1)

})


test_that("gt_sample yields expected mean GT", {
  set.seed(1)
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  gt <- gt_sample(p, n = 10000)
  ## large tolerance because of discretisation
  expect_true(abs(mean(gt) - 6.3) < 0.5) # value from Bi et al.
})


test_that("gt_sample yields expected invalid input errors", {
  p_base <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  p <- p_base
  p$I_C_2_transmission <- 0.1
  expect_error(
    gt_sample(p, n = 1000),
    "gt_sample does not allow transmission from I_C_2",
    fixed = TRUE)
  p <- p_base
  p$I_P_transmission <- 0.1
  expect_error(
    gt_sample(p, n = 1000),
    "gt_sample does not allow
    I_P_transmission !=1 or I_C_1_transmission != 1",
    fixed = TRUE)
  p <- p_base
  p$I_C_1_transmission <- 0.1
  expect_error(
    gt_sample(p, n = 1000),
    "gt_sample does not allow
    I_P_transmission !=1 or I_C_1_transmission != 1",
    fixed = TRUE)
  p <- p_base
  p$k_A <- 2
  expect_error(
    gt_sample(p, n = 1000),
    "gt_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1",
    fixed = TRUE)
  p <- p_base
  p$k_P <- 2
  expect_error(
    gt_sample(p, n = 1000),
    "gt_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1",
    fixed = TRUE)
  p <- p_base
  p$k_C_1 <- 2
  expect_error(
    gt_sample(p, n = 1000),
    "gt_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1",
    fixed = TRUE)

})

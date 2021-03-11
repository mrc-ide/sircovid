test_that("Can calculate EpiEstim Rt", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)
  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)
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

  rt_EpiEstim <- carehomes_EpiEstim_Rt_trajectories(steps, inc, p,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  #### General patterns
  ## dimension of Rt should be 3 rows and length(steps) - 1 cols
  ## the -1 because EpiEstim only start estimation at 2nd time step
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(length(steps) - 1, 4))
  ## Rt at the end should be < 1
  expect_true(all(rt_EpiEstim$Rt_summary[nrow(rt_EpiEstim$Rt_summary)] < 1))
  ## because of the priors with mean 1 we would expect EpiEstim Rt
  ## at the end to be higher than eff_Rt_all
  expect_true(
    last(rt$eff_Rt_all) <
      rt_EpiEstim$Rt_summary[nrow(rt_EpiEstim$Rt_summary), "mean_R"])
  ## Rt at the start should be > 1
  ## (but there will be uncertainty so looking at the mean)
  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary[, "mean_R"])))
  expect_true(rt_EpiEstim$Rt_summary[first_non_NA_idx, "mean_R"] > 1)

  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[[first_non_NA_idx, "mean_R"]], 4.8,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[[first_non_NA_idx, "2.5%"]], 1.5,
               tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[[nrow(rt_EpiEstim$Rt_summary), "mean_R"]], 0.1,
    tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[[nrow(rt_EpiEstim$Rt_summary), "97.5%"]], 0.1,
    tolerance = .1)

})


test_that("Can calculate EpiEstim Rt when no transmission in carehomes", {
  skip_if_not_installed("EpiEstim")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            beta_value = 0.15)
  mean_N_tot <- mean(p$N_tot)
  mean_m <- mean(p$m)
  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)
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

  rt_EpiEstim <- carehomes_EpiEstim_Rt_trajectories(steps, inc, p,
                                                    sliding_window_ndays = 1)

  rt <- carehomes_Rt_trajectories(steps, S, p)

  #### General patterns
  ## dimension of Rt should be 3 rows and length(steps) - 1 cols
  ## the -1 because EpiEstim only start estimation at 2nd time step
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(length(steps) - 1, 4))
  ## Rt at the end should be < 1
  expect_true(all(rt_EpiEstim$Rt_summary[nrow(rt_EpiEstim$Rt_summary), ] < 1))
  ## because of the priors with mean 1 we would expect EpiEstim Rt
  ## at the end to be higher than eff_Rt_all
  expect_true(
    last(rt$eff_Rt_all) <
      rt_EpiEstim$Rt_summary[nrow(rt_EpiEstim$Rt_summary), "mean_R"])
  ## Rt at the start should be > 1
  ## (but there will be uncertainty so looking at the mean)
  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary[, "mean_R"])))
  expect_true(rt_EpiEstim$Rt_summary[first_non_NA_idx, "mean_R"] > 1)
  ## Rt at the start should be similar between methods:
  expect_true(
    rt$eff_Rt_all[1] > rt_EpiEstim$Rt_summary[first_non_NA_idx, "2.5%"])
  expect_true(
    rt$eff_Rt_all[1] < rt_EpiEstim$Rt_summary[first_non_NA_idx, "97.5%"])

  ## Rt stays constant for a while but precision improves for EpiEstim
  expect_true(rt$eff_Rt_all[1] > rt_EpiEstim$Rt_summary[20, "2.5%"])
  expect_true(rt$eff_Rt_all[1] < rt_EpiEstim$Rt_summary[20, "97.5%"])
  expect_true(abs(rt$eff_Rt_all[1] - rt_EpiEstim$Rt_summary[20, "50%"]) < 0.1)

  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[[35, "mean_R"]], 6.5,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[[35, "2.5%"]], 5.8,
               tolerance = .1)
  expect_equal(rt$eff_Rt_all[35], 6.5,
               tolerance = .1)

  expect_equal(rt_EpiEstim$Rt_summary[[45, "mean_R"]], 0.4,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[[45, "2.5%"]], 0.1,
               tolerance = .1)
  expect_equal(rt$eff_Rt_all[45], 0.6,
               tolerance = .1)

})


test_that("gt_sample yields expected mean GT", {
  set.seed(1)
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)
  gt <- gt_sample(p, n = 10000)
  ## large tolerance because of discretisation
  expect_true(abs(mean(gt) - 6.3) < 0.5) # value from Bi et al.
})


test_that("gt_sample yields expected invalid input errors", {
  p_base <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  expect_error(
    gt_sample(p_base, n = 1000),
    "gt_sample does not allow p_C to vary by age",
    fixed = TRUE)

  ## Fix p_C across age groups for the rest of the test
  p_base$p_C <- rep(0.6, 19)

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

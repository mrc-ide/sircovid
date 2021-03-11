context("carehomes (Rt)")

test_that("Can calculate Rt", {
  ## This does not really require that we have run the real model, as
  ## everything is deterministic. So here we'll generate some data
  ## from the model and use it so that we can compare against some
  ## exact numbers.
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  res <- carehomes_Rt(steps, y[, 1, ], p)
  expect_equal(res, d$outputs$rt_1)
  res_all <- carehomes_Rt_trajectories(steps, y, p)
  expect_equal(res_all, d$outputs$rt_all)

  ## Beta is returned, results correct length
  expect_identical(res$beta, rep(p$beta_step, length(steps)))
  expect_true(all(lengths(res) == length(steps)))

  ## Check the Rt calculation (from eff_Rt)
  expect_true(diff(range(res$Rt_all)) < 1e-7)
  expect_true(diff(range(res$Rt_general)) < 1e-7)

  ## Effective Rt lower than Rt
  expect_true(all(res$Rt_all - res$eff_Rt_all > -1e-7))
  expect_true(all(res$Rt_general - res$eff_Rt_general > -1e-7))

  ## General population effective Rt lower than total
  expect_true(all(res$eff_Rt_all >= res$eff_Rt_general))

  expect_equal(names(res), names(res_all))
  for (nm in names(res)) {
    expect_equal(res[[nm]], res_all[[nm]][, 1, drop = TRUE])
  }

  ## Date is returned
  expect_equal(res$date, res$step * p$dt)
})


test_that("validate inputs in Rt calculation", {
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_Rt(steps, y[-1, 1, ], p),
    "Expected 'S' to have 19 rows = 19 groups x 1 vaccine classes",
    fixed = TRUE)
  expect_error(
    carehomes_Rt(steps, y[, 1, -1], p),
    "Expected 'S' to have 85 columns, following 'step'",
    fixed = TRUE)
})


test_that("validate inputs in Rt trajectories calculation", {
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_Rt_trajectories(steps, y[, 1, ], p),
    "Expected a 3d array of 'S'",
    fixed = TRUE)
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , -1], p),
    "Expected 3rd dimension of 'S' to have length 85, following 'step'",
    fixed = TRUE)
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , ], list(p)),
    "Expected 2nd dimension of 'S' to have length 1, following 'pars'")
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , ], p, shared_parameters = FALSE),
    "If not using shared parameters, expected a unnamed list for 'pars'")
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , ], list(p),
                              shared_parameters = TRUE),
    "If using shared parameters, expected a named list for 'pars'")
})


test_that("Can set initial time", {
  ## This test also checks that we can alter parameter inputs
  d <- reference_data_rt()

  steps <- d$inputs$steps
  y <- d$inputs$y

  step0 <- seq(steps[[1]], by = 1, length.out = ncol(y))

  p <- rep(list(d$inputs$p), ncol(y))
  for (i in seq_along(p)) {
    p[[i]]$initial_step <- step0[[i]]
  }

  res1 <- carehomes_Rt_trajectories(steps, y, p,
                                    initial_step_from_parameters = FALSE)
  expect_equal(res1$step, matrix(steps, length(steps), ncol(y)))

  res2 <- carehomes_Rt_trajectories(steps, y, p,
                                    initial_step_from_parameters = TRUE)
  expect_equal(res2$step[1, ], step0)
  expect_equal(res2$step[-1, ], res1$step[-1, ])
})


test_that("Can vary beta over time", {
  d <- reference_data_rt()

  steps <- d$inputs$steps
  y <- d$inputs$y

  p <- d$inputs$p
  dt <- p$dt
  initial_date <- p$initial_step * dt
  beta_date <- initial_date + c(0, 21, 62)
  beta_value <- p$beta_step * c(1, 0.5, 0.8)
  p$beta_step <- sircovid_parameters_beta(beta_date, beta_value, dt)
  p <- rep(list(p), ncol(y))
  p[[3]]$beta_step <- p[[3]]$beta_step / 2

  res <- carehomes_Rt_trajectories(steps, y, p,
                                   initial_step_from_parameters = FALSE)

  expect_equal(
    res$beta,
    cbind(sircovid_parameters_beta_expand(steps, p[[1]]$beta_step),
          sircovid_parameters_beta_expand(steps, p[[2]]$beta_step),
          sircovid_parameters_beta_expand(steps, p[[3]]$beta_step)))

  ## Check the Rt calculation (from eff_Rt) - compare the first test
  expect_true(length(unique(res$Rt_all)) > 1)
  expect_true(length(unique(res$Rt_general)) > 1)
  expect_true(all(res$Rt_all - res$eff_Rt_all > -1e-7))
  expect_true(all(res$Rt_general >= res$eff_Rt_general))
})


test_that("can filter Rt to wanted types", {
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_mapequal(
    carehomes_Rt(steps, y[, 1, ], p, type = "eff_Rt_general"),
    d$outputs$rt_1[c("step", "date", "beta", "eff_Rt_general")])
  expect_mapequal(
    carehomes_Rt(steps, y[, 1, ], p, type = c("eff_Rt_general", "Rt_general")),
    d$outputs$rt_1[c("step", "date", "beta", "eff_Rt_general", "Rt_general")])

  expect_mapequal(
    carehomes_Rt_trajectories(steps, y, p, type = "eff_Rt_all"),
    d$outputs$rt_all[c("step", "date", "beta", "eff_Rt_all")])
  expect_mapequal(
    carehomes_Rt_trajectories(steps, y, p, type = c("eff_Rt_all", "Rt_all")),
    d$outputs$rt_all[c("step", "date", "beta", "eff_Rt_all", "Rt_all")])
})


test_that("can't compute Rt for unknown types", {
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_Rt(steps, y[, 1, ], p, type = "max_Rt_general"),
    "Unknown R type 'max_Rt_general', must match '")
  expect_error(
    carehomes_Rt_trajectories(steps, y, p, type = "max_Rt_general"),
    "Unknown R type 'max_Rt_general', must match '")
  expect_error(
    carehomes_Rt(steps, y[, 1, ], p, type = c("eff_Rt_general", "rt_general")),
    "Unknown R type 'rt_general', must match '")
})

test_that("Can interpolate Rt with step changes", {
  dat <- reference_data_mcmc()
  rt <- local({
    p <- lapply(seq_len(nrow(dat$pars)), function(i)
      dat$predict$transform(dat$pars[i, ]))
    i <- grep("S_", rownames(dat$trajectories$state))
    S <- dat$trajectories$state[i, , ]
    carehomes_Rt_trajectories(dat$trajectories$step, S, p)
  })

  future <- list(
    "2020-04-15" = future_Rt(1.5, "2020-03-10"),
    "2020-05-01" = future_Rt(0.5, "2020-03-10"),
    "2020-05-15" = future_Rt(2, "2020-03-10"))
  res <- future_relative_beta(future, rt$date[, 1], rt$Rt_general)
  baseline <- add_future_betas(dat, rt, future)

  events <- sircovid_simulate_events("2020-03-30", "2020-06-01", NULL)
  p <- lapply(seq_len(nrow(baseline$pars)), function(i)
    baseline$predict$transform(baseline$pars[i, ]))
  ans <- sircovid_simulate(carehomes, baseline$state, p, events,
                           index = dat$predict$index)

  ## Work out our critical dates so that we can start interpolation:
  step <- attr(ans, "step")
  S <- ans[grep("S_", rownames(ans)), , ]

  crit_dates <- sircovid_date(names(future))

  set.seed(1)
  rt_cmp <- carehomes_Rt_trajectories(step, S, p,
                                      initial_step_from_parameters = FALSE)

  ## Only interpolate if "every" is given:
  set.seed(1)
  expect_identical(
    carehomes_Rt_trajectories(step, S, p,
                              initial_step_from_parameters = FALSE,
                              interpolate_min = 3),
    rt_cmp)


  ## Then compute the Rt values with interpolation
  rt_int_2 <- carehomes_Rt_trajectories(step, S, p,
                                        initial_step_from_parameters = FALSE,
                                        interpolate_every = 2,
                                        interpolate_min = 3,
                                        interpolate_critical_dates = crit_dates)
  rt_int_7 <- carehomes_Rt_trajectories(step, S, p,
                                        initial_step_from_parameters = FALSE,
                                        interpolate_every = 7,
                                        interpolate_min = 3,
                                        interpolate_critical_dates = crit_dates)
  rt_int_14 <- carehomes_Rt_trajectories(step, S, p,
                                         initial_step_from_parameters = FALSE,
                                         interpolate_every = 14,
                                         interpolate_min = 1,
                                         interpolate_critical_dates =
                                           crit_dates)
  ## check the error is small
  tol <- 0.05
  # for interpolation every 2 days
  expect_true(all(abs(rt_cmp$eff_Rt_all - rt_int_2$eff_Rt_all) < tol))
  expect_true(all(abs(rt_cmp$eff_Rt_general - rt_int_2$eff_Rt_general) < tol))
  expect_true(all(abs(rt_cmp$Rt_all - rt_int_2$Rt_all) < tol))
  expect_true(all(abs(rt_cmp$Rt_general - rt_int_2$Rt_general) < tol))
  # for interpolation every 7 days
  expect_true(all(abs(rt_cmp$eff_Rt_all - rt_int_7$eff_Rt_all) < tol))
  expect_true(all(abs(rt_cmp$eff_Rt_general - rt_int_7$eff_Rt_general) < tol))
  expect_true(all(abs(rt_cmp$Rt_all - rt_int_7$Rt_all) < tol))
  expect_true(all(abs(rt_cmp$Rt_general - rt_int_7$Rt_general) < tol))
  # have to increase tolerance dramatically for every 14 days
  tol2 <- 0.5
  expect_true(all(abs(rt_cmp$eff_Rt_all - rt_int_14$eff_Rt_all) < tol2))
  expect_true(all(abs(rt_cmp$eff_Rt_general - rt_int_14$eff_Rt_general) < tol2))
  expect_true(all(abs(rt_cmp$Rt_all - rt_int_14$Rt_all) < tol2))
  expect_true(all(abs(rt_cmp$Rt_general - rt_int_14$Rt_general) < tol2))
})


test_that("Parameters affect Rt as expected", {

  ## Note that m_CHW and m_CHR have been changed from defaults to avoid
  ## having all care home residents infected
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            m_CHW = 3e-6, m_CHR = 3e-6)

  ## set the following parameters to non-zero values to allow related parameters
  ## to have an effect on Rt
  p$hosp_transmission <- 0.05
  p$ICU_transmission <- 0.05
  p$G_D_transmission <- 0.05
  p$I_C_2_transmission <- 0.5

  ## allow some deaths in the community so that G_D parameters affect Rt_general
  p$psi_G_D[15:17] <- 0.5

  np <- 1L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  helper <- function(par_name, par_value_lower_Rt, par_value_higher_Rt) {

    p[[par_name]] <- par_value_lower_Rt
    rt_lower <- carehomes_Rt(steps, y[, 1, ], p)

    p[[par_name]] <- par_value_higher_Rt
    rt_higher <- carehomes_Rt(steps, y[, 1, ], p)

    expect_true(all(rt_lower$Rt_all < rt_higher$Rt_all))
    expect_true(all(rt_lower$eff_Rt_all < rt_higher$eff_Rt_all))
    expect_true(all(rt_lower$Rt_general < rt_higher$Rt_general))
    expect_true(all(rt_lower$eff_Rt_general < rt_higher$eff_Rt_general))
  }

  helper("I_A_transmission", 0, 1)
  helper("I_P_transmission", 0, 1)
  helper("I_C_1_transmission", 0, 1)
  helper("I_C_2_transmission", 0, 1)
  helper("hosp_transmission", 0, 1)
  helper("ICU_transmission", 0, 1)
  helper("G_D_transmission", 0, 1)

  helper("gamma_A", Inf, 1)
  helper("gamma_P", Inf, 1)
  helper("gamma_C_1", Inf, 1)
  helper("gamma_C_2", Inf, 1)
  helper("gamma_H_D", Inf, 1)
  helper("gamma_H_R", Inf, 1)
  helper("gamma_ICU_pre", Inf, 1)
  helper("gamma_ICU_D", Inf, 1)
  helper("gamma_ICU_W_D", Inf, 1)
  helper("gamma_ICU_W_R", Inf, 1)
  helper("gamma_G_D", Inf, 1)

  helper("k_A", 1, 2)
  helper("k_P", 1, 2)
  helper("k_C_1", 1, 2)
  helper("k_C_2", 1, 2)
  helper("k_H_D", 1, 2)
  helper("k_H_R", 1, 2)
  helper("k_ICU_pre", 1, 2)
  helper("k_ICU_D", 1, 2)
  helper("k_ICU_W_D", 1, 2)
  helper("k_ICU_W_R", 1, 2)
  helper("k_G_D", 1, 2)

})


test_that("Can calculate EpiEstim Rt", {
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
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(4, length(steps) - 1))
  ## Rt at the end should be < 1
  expect_true(all(rt_EpiEstim$Rt_summary[, ncol(rt_EpiEstim$Rt_summary)] < 1))
  ## because of the priors with mean 1 we would expect EpiEstim Rt
  ## at the end to be higher than eff_Rt_all
  expect_true(
    last(rt$eff_Rt_all) <
      rt_EpiEstim$Rt_summary["mean_R", ncol(rt_EpiEstim$Rt_summary)])
  ## Rt at the start should be > 1
  ## (but there will be uncertainty so looking at the mean)
  first_non_NA_idx <- min(which(!is.na(rt_EpiEstim$Rt_summary["mean_R", ])))
  expect_true(rt_EpiEstim$Rt_summary["mean_R", first_non_NA_idx] > 1)

  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", first_non_NA_idx]], 4.8,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", first_non_NA_idx]], 1.5,
               tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[["mean_R", ncol(rt_EpiEstim$Rt_summary)]], 0.1,
               tolerance = .1)
  expect_equal(
    rt_EpiEstim$Rt_summary[["97.5%", ncol(rt_EpiEstim$Rt_summary)]], 0.1,
               tolerance = .1)

})


test_that("Can calculate EpiEstim Rt when no transmission in carehomes", {
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
  expect_equal(dim(rt_EpiEstim$Rt_summary), c(4, length(steps) - 1))
  ## Rt at the end should be < 1
  expect_true(all(rt_EpiEstim$Rt_summary[, ncol(rt_EpiEstim$Rt_summary)] < 1))
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
  expect_true(abs(rt$eff_Rt_all[1] - rt_EpiEstim$Rt_summary["50%", 20]) < 0.1)

  #### Check a few values
  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", 35]], 6.5,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", 35]], 5.8,
               tolerance = .1)
  expect_equal(rt$eff_Rt_all[35], 6.5,
               tolerance = .1)

  expect_equal(rt_EpiEstim$Rt_summary[["mean_R", 45]], 0.4,
               tolerance = .1)
  expect_equal(rt_EpiEstim$Rt_summary[["2.5%", 45]], 0.1,
               tolerance = .1)
  expect_equal(rt$eff_Rt_all[45], 0.6,
               tolerance = .1)

})


test_that("draw_one_GT_sample yields expected mean GT", {
  set.seed(1)
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  ## Fix p_C across age groups for the rest of the test
  p$p_C <- rep(0.6, 19)
  gt <- draw_one_GT_sample(p, n = 10000)
  ## large tolerance because of discretisation
  expect_true(abs(mean(gt) - 6.3) < 0.5) # value from Bi et al.
})


test_that("draw_one_GT_sample yields expected invalid input errors", {
  p_base <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  expect_error(
    draw_one_GT_sample(p_base, n = 1000),
    "draw_one_GT_sample does not allow p_C to vary by age",
    fixed = TRUE)

  ## Fix p_C across age groups for the rest of the test
  p_base$p_C <- rep(0.6, 19)

  p <- p_base
  p$I_C_2_transmission <- 0.1
  expect_error(
    draw_one_GT_sample(p, n = 1000),
    "draw_one_GT_sample does not allow transmission from I_C_2",
    fixed = TRUE)
  p <- p_base
  p$I_P_transmission <- 0.1
  expect_error(
    draw_one_GT_sample(p, n = 1000),
    "draw_one_GT_sample does not allow
    I_P_transmission !=1 or I_C_1_transmission != 1",
    fixed = TRUE)
  p <- p_base
  p$I_C_1_transmission <- 0.1
  expect_error(
    draw_one_GT_sample(p, n = 1000),
    "draw_one_GT_sample does not allow
    I_P_transmission !=1 or I_C_1_transmission != 1",
    fixed = TRUE)
  p <- p_base
  p$k_A <- 2
  expect_error(
    draw_one_GT_sample(p, n = 1000),
    "draw_one_GT_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1",
    fixed = TRUE)
  p <- p_base
  p$k_P <- 2
  expect_error(
    draw_one_GT_sample(p, n = 1000),
    "draw_one_GT_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1",
    fixed = TRUE)
  p <- p_base
  p$k_C_1 <- 2
  expect_error(
    draw_one_GT_sample(p, n = 1000),
    "draw_one_GT_sample does not allow k_A > 1, k_P > 1 or k_C_1 > 1",
    fixed = TRUE)

})

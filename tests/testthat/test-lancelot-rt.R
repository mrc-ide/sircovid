context("lancelot (Rt)")

test_that("Can calculate Rt", {
  ## This does not really require that we have run the real model, as
  ## everything is deterministic. So here we'll generate some data
  ## from the model and use it so that we can compare against some
  ## exact numbers.
  d <- reference_data_lancelot_rt()

  p <- d$inputs$p
  time <- d$inputs$time
  y <- d$inputs$y

  res <- lancelot_Rt(time, y[, 1, ], p)
  expect_equal(unclass(res), unclass(d$outputs$rt_1))
  res_all <- lancelot_Rt_trajectories(time, y, p)
  expect_equal(unclass(res_all), unclass(d$outputs$rt_all))

  ## Beta is returned, results correct length
  expect_identical(res$beta, rep(p$beta_step, length(time)))
  expect_vector_equal(lengths(res), length(time))

  ## Check the Rt calculation (from eff_Rt)
  expect_true(diff(range(res$Rt_all)) < 1e-7)
  expect_true(diff(range(res$Rt_general)) < 1e-7)

  ## Effective Rt lower than Rt
  expect_vector_gt(res$Rt_all, res$eff_Rt_all, tol = 1e-7)
  expect_vector_gt(res$Rt_general, res$eff_Rt_general, tol = 1e-7)

  ## General population effective Rt lower than total
  expect_vector_gte(res$eff_Rt_all, res$eff_Rt_general)

  expect_equal(names(res), names(res_all))
  for (nm in names(res)) {
    expect_equal(drop(res[[nm]]), res_all[[nm]][, 1, drop = TRUE])
  }

  ## Date is returned
  expect_equal(res$date, res$time * p$dt)
})


test_that("validate inputs in Rt calculation", {
  d <- reference_data_lancelot_rt()

  p <- d$inputs$p
  time <- d$inputs$time
  y <- d$inputs$y

  expect_error(
    lancelot_Rt(time, y[-1, 1, ], p),
    "Expected 'S' to have 19 rows = 19 groups x 1 vaccine classes",
    fixed = TRUE)
  expect_error(
    lancelot_Rt(time, y[, 1, -1], p),
    "Expected 'S' to have 123 columns, following 'time'",
    fixed = TRUE)
})


test_that("validate inputs in Rt trajectories calculation", {
  d <- reference_data_lancelot_rt()

  p <- d$inputs$p
  time <- d$inputs$time
  y <- d$inputs$y

  expect_error(
    lancelot_Rt_trajectories(time, y[, 1, ], p),
    "Expected a 3d array of 'S'",
    fixed = TRUE)
  expect_error(
    lancelot_Rt_trajectories(time, y[, , -1], p),
    "Expected 3rd dimension of 'S' to have length 123, following 'time'",
    fixed = TRUE)
  expect_error(
    lancelot_Rt_trajectories(time, y[, , ], list(p)),
    "Expected 2nd dimension of 'S' to have length 1, following 'pars'")
  expect_error(
    lancelot_Rt_trajectories(time, y[, , ], p, shared_parameters = FALSE),
    "If not using shared parameters, expected a unnamed list for 'pars'")
  expect_error(
    lancelot_Rt_trajectories(time, y[, , ], list(p),
                             shared_parameters = TRUE),
    "If using shared parameters, expected a named list for 'pars'")
})


test_that("Can set initial time", {
  ## This test also checks that we can alter parameter inputs
  d <- reference_data_lancelot_rt()

  time <- d$inputs$time
  y <- d$inputs$y

  time0 <- seq(time[[1]], by = 1, length.out = ncol(y))

  p <- rep(list(d$inputs$p), ncol(y))
  for (i in seq_along(p)) {
    p[[i]]$initial_time <- time0[[i]]
  }

  res1 <- lancelot_Rt_trajectories(time, y, p,
                                   initial_time_from_parameters = FALSE)
  expect_equal(res1$time, matrix(time, length(time), ncol(y)))

  res2 <- lancelot_Rt_trajectories(time, y, p,
                                   initial_time_from_parameters = TRUE)
  expect_equal(res2$time[1, ], time0)
  expect_equal(res2$time[-1, ], res1$time[-1, ])
})


test_that("Can vary beta over time", {
  d <- reference_data_lancelot_rt()

  time <- d$inputs$time
  y <- d$inputs$y

  p <- d$inputs$p
  dt <- p$dt
  initial_date <- 0
  beta_date <- initial_date + c(0, 21, 62)
  beta_value <- p$beta_step * c(1, 0.5, 0.8)
  p$beta_step <- sircovid_parameters_piecewise_linear(beta_date,
                                                      beta_value, dt)
  p <- rep(list(p), ncol(y))
  p[[3]]$beta_step <- p[[3]]$beta_step / 2

  res <- lancelot_Rt_trajectories(time, y, p,
                                  initial_time_from_parameters = FALSE)

  expect_equal(
    res$beta,
    cbind(sircovid_parameters_expand_step(time, p[[1]]$beta_step),
          sircovid_parameters_expand_step(time, p[[2]]$beta_step),
          sircovid_parameters_expand_step(time, p[[3]]$beta_step)))

  ## Check the Rt calculation (from eff_Rt) - compare the first test
  expect_true(length(unique(res$Rt_all)) > 1)
  expect_true(length(unique(res$Rt_general)) > 1)
  expect_vector_gt(res$Rt_all, res$eff_Rt_all, tol = 1e-7)
  expect_vector_gte(res$Rt_general, res$eff_Rt_general, tol = 1e-7)
})


test_that("can filter Rt to wanted types", {
  d <- reference_data_lancelot_rt()

  p <- d$inputs$p
  time <- d$inputs$time
  y <- d$inputs$y

  expect_mapequal(
    lancelot_Rt(time, y[, 1, ], p, type = "eff_Rt_general"),
    d$outputs$rt_1[c("time", "date", "beta", "eff_Rt_general")])
  expect_mapequal(
    lancelot_Rt(time, y[, 1, ], p, type = c("eff_Rt_general", "Rt_general")),
    d$outputs$rt_1[c("time", "date", "beta", "eff_Rt_general", "Rt_general")])

  expect_mapequal(
    lancelot_Rt_trajectories(time, y, p, type = "eff_Rt_all"),
    d$outputs$rt_all[c("time", "date", "beta", "eff_Rt_all")])
  expect_mapequal(
    lancelot_Rt_trajectories(time, y, p, type = c("eff_Rt_all", "Rt_all")),
    d$outputs$rt_all[c("time", "date", "beta", "eff_Rt_all", "Rt_all")])
})


test_that("can't compute Rt for unknown types", {
  d <- reference_data_lancelot_rt()

  p <- d$inputs$p
  time <- d$inputs$time
  y <- d$inputs$y

  expect_error(
    lancelot_Rt(time, y[, 1, ], p, type = "max_Rt_general"),
    "Unknown R type 'max_Rt_general', must match '")
  expect_error(
    lancelot_Rt_trajectories(time, y, p, type = "max_Rt_general"),
    "Unknown R type 'max_Rt_general', must match '")
  expect_error(
    lancelot_Rt(time, y[, 1, ], p, type = c("eff_Rt_general", "rt_general")),
    "Unknown R type 'rt_general', must match '")
})


test_that("Can interpolate Rt with time changes", {
  skip("Needs updating without add_future_betas/sircovid_simulate*")
  dat <- reference_data_lancelot_mcmc()
  rt <- local({
    p <- lapply(seq_len(nrow(dat$pars)), function(i)
      dat$predict$transform(dat$pars[i, ]))
    i <- grep("S_", rownames(dat$trajectories$state))
    S <- dat$trajectories$state[i, , ]
    lancelot_Rt_trajectories(dat$trajectories$time, S, p)
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
  ans <- sircovid_simulate(lancelot, baseline$state, p, events,
                           index = dat$predict$index)

  ## Work out our critical dates so that we can start interpolation:
  time <- attr(ans, "time")
  S <- ans[grep("S_", rownames(ans)), , ]

  crit_dates <- sircovid_date(names(future))

  set.seed(1)
  rt_cmp <- lancelot_Rt_trajectories(time, S, p,
                                     initial_time_from_parameters = FALSE)

  ## Only interpolate if "every" is given:
  set.seed(1)
  expect_identical(
    lancelot_Rt_trajectories(time, S, p,
                             initial_time_from_parameters = FALSE,
                             interpolate_min = 3),
    rt_cmp)


  ## Then compute the Rt values with interpolation
  rt_int_2 <- lancelot_Rt_trajectories(time, S, p,
                                       initial_time_from_parameters = FALSE,
                                       interpolate_every = 2,
                                       interpolate_min = 3,
                                       interpolate_critical_dates = crit_dates)
  rt_int_7 <- lancelot_Rt_trajectories(time, S, p,
                                       initial_time_from_parameters = FALSE,
                                       interpolate_every = 7,
                                       interpolate_min = 3,
                                       interpolate_critical_dates = crit_dates)
  rt_int_14 <- lancelot_Rt_trajectories(time, S, p,
                                        initial_time_from_parameters = FALSE,
                                        interpolate_every = 14,
                                        interpolate_min = 1,
                                        interpolate_critical_dates =
                                          crit_dates)
  ## check the error is small
  tol <- 0.05
  # for interpolation every 2 days
  expect_vector_equal(rt_cmp$eff_Rt_all, rt_int_2$eff_Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$eff_Rt_general, rt_int_2$eff_Rt_general,
                      tol = tol)
  expect_vector_equal(rt_cmp$Rt_all, rt_int_2$Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$Rt_general, rt_int_2$Rt_general, tol = tol)
  # for interpolation every 7 days
  expect_vector_equal(rt_cmp$eff_Rt_all, rt_int_7$eff_Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$eff_Rt_general,
                      rt_int_7$eff_Rt_general, tol = tol)
  expect_vector_equal(rt_cmp$Rt_all, rt_int_7$Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$Rt_general, rt_int_7$Rt_general, tol = tol)
  # have to increase tolerance dramatically for every 14 days
  tol2 <- 0.5
  expect_vector_equal(rt_cmp$eff_Rt_all, rt_int_14$eff_Rt_all, tol = tol2)
  expect_vector_equal(rt_cmp$eff_Rt_general,
                      rt_int_14$eff_Rt_general, tol = tol2)
  expect_vector_equal(rt_cmp$Rt_all, rt_int_14$Rt_all, tol = tol2)
  expect_vector_equal(rt_cmp$Rt_general, rt_int_14$Rt_general, tol = tol2)
})


test_that("Parameters affect Rt as expected", {

  ## Note that m_CHW and m_CHR have been changed from defaults to avoid
  ## having all care home residents infected
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           m_CHW = 3e-6, m_CHR = 3e-6)

  ## set the following parameters to non-zero values to allow related parameters
  ## to have an effect on Rt
  p$I_C_2_transmission <- 0.5

  np <- 1L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  time <- seq(0, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(time)

  helper <- function(par_name, par_value_lower_Rt, par_value_higher_Rt) {

    p[[par_name]] <- par_value_lower_Rt
    rt_lower <- lancelot_Rt(time, y[, 1, ], p)

    p[[par_name]] <- par_value_higher_Rt
    rt_higher <- lancelot_Rt(time, y[, 1, ], p)

    expect_vector_lt(rt_lower$Rt_all, rt_higher$Rt_all)
    expect_vector_lt(rt_lower$eff_Rt_all, rt_higher$eff_Rt_all)
    expect_vector_lt(rt_lower$Rt_general, rt_higher$Rt_general)
    expect_vector_lt(rt_lower$eff_Rt_general, rt_higher$eff_Rt_general)
  }

  helper("I_A_transmission", 0, 1)
  helper("I_P_transmission", 0, 1)
  helper("I_C_1_transmission", 0, 1)
  helper("I_C_2_transmission", 0, 1)

  helper("gamma_A_step", Inf, 1)
  helper("gamma_P_step", Inf, 1)
  helper("gamma_C_1_step", Inf, 1)
  helper("gamma_C_2_step", Inf, 1)

  helper("k_A", 1, 2)
  helper("k_P", 1, 2)
  helper("k_C_1", 1, 2)
  helper("k_C_2", 1, 2)
})

test_that("p_C capped correctly at 1 in Rt calculation", {
  ## p_C is the only severity probability that impacts the Rt calculation

  ## Note that m_CHW and m_CHR have been changed from defaults to avoid
  ## having all care home residents infected
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           m_CHW = 3e-6, m_CHR = 3e-6)

  np <- 1L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  time <- seq(0, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(time)

  p$p_C_step[, ] <- 1
  rt1 <- lancelot_Rt(time, y[, 1, ], p)

  p$p_C_step[, ] <- 1.5
  rt2 <- lancelot_Rt(time, y[, 1, ], p)

  expect_equal(rt1$Rt_all, rt2$Rt_all)
  expect_equal(rt1$eff_Rt_all, rt2$eff_Rt_all)
  expect_equal(rt1$Rt_general, rt2$Rt_general)
  expect_equal(rt1$eff_Rt_general, rt2$eff_Rt_general)

})

test_that("Cannot calculate Rt for non-zero transmissions", {
  d <- reference_data_lancelot_rt()
  time <- d$inputs$time
  y <- d$inputs$y

  p <- d$inputs$p
  p$hosp_transmission <- 1
  expect_error(lancelot_Rt(time, y[, 1, ], p), "Cannot currently compute")

  p <- d$inputs$p
  p$G_D_transmission <- 1
  expect_error(lancelot_Rt(time, y[, 1, ], p), "Cannot currently compute")

  p <- d$inputs$p
  p$ICU_transmission <- 1
  expect_error(lancelot_Rt(time, y[, 1, ], p), "Cannot currently compute")
})

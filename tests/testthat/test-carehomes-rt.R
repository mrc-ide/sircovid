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
  expect_true(length(unique(res$Rt_all)) == 1)
  expect_true(length(unique(res$Rt_general)) == 1)

  ## Effective Rt lower than Rt
  expect_true(all(res$Rt_all >= res$eff_Rt_all))
  expect_true(all(res$Rt_general >= res$eff_Rt_general))

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
  expect_true(all(res$Rt_all >= res$eff_Rt_all))
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

  rt_cmp <- carehomes_Rt_trajectories(step, S, p,
                                      initial_step_from_parameters = FALSE)

  ## Only interpolate if "every" is given:
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
                                        interpolate_critical_dates = crit_dates)
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

context("carehomes (IFR_t)")

test_that("Can calculate IFR_t", {
  ## This does not really require that we have run the real model, as
  ## everything is deterministic. So here we'll generate some data
  ## from the model and use it so that we can compare against some
  ## exact numbers.
  d <- reference_data_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  res <- carehomes_ifr_t(steps, y[, 1, ], p)
  expect_equal(res, d$outputs$ifr_t_1)
  res_all <- carehomes_ifr_t_trajectories(steps, y, p)
  expect_equal(res_all, d$outputs$ifr_t_all)

  ## Results correct length
  expect_true(all(lengths(res) == length(steps)))

  expect_equal(names(res), names(res_all))
  for (nm in names(res)) {
    expect_equal(res[[nm]], res_all[[nm]][, 1, drop = TRUE])
  }

  ## Date is returned
  expect_equal(res$date, res$step * p$dt)
})


test_that("validate inputs in ifr_t calculation", {
  d <- reference_data_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_ifr_t(steps, y[-1, 1, ], p),
    "Expected 'infections_inc' to have 19 rows = 19 groups x 1 vacc classes",
    fixed = TRUE)
  expect_error(
    carehomes_ifr_t(steps, y[, 1, -1], p),
    "Expected 'infections_inc' to have 85 cols, following 'step'",
    fixed = TRUE)
})


test_that("validate inputs in rt trajectories calculation", {
  d <- reference_data_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_ifr_t_trajectories(steps, y[, 1, ], p),
    "Expected a 3d array of 'infections_inc'",
    fixed = TRUE)
  expect_error(
    carehomes_ifr_t_trajectories(steps, y[, , -1], p),
    "Expected 3rd dim of 'infections_inc' to have length 85, given 'step'",
    fixed = TRUE)
  expect_error(
    carehomes_ifr_t_trajectories(steps, y[, , ], list(p)),
    "Expected 2nd dim of 'infections_inc' to have length 1, given 'pars'")
  expect_error(
    carehomes_ifr_t_trajectories(steps, y[, , ], p, shared_parameters = FALSE),
    "If not using shared parameters, expected a unnamed list for 'pars'")
  expect_error(
    carehomes_ifr_t_trajectories(steps, y[, , ], list(p),
                              shared_parameters = TRUE),
    "If using shared parameters, expected a named list for 'pars'")
})


test_that("Can set initial time", {
  ## This test also checks that we can alter parameter inputs
  d <- reference_data_ifr_t()

  steps <- d$inputs$steps
  y <- d$inputs$y

  step0 <- seq(steps[[1]], by = 1, length.out = ncol(y))

  p <- rep(list(d$inputs$p), ncol(y))
  for (i in seq_along(p)) {
    p[[i]]$initial_step <- step0[[i]]
  }

  res1 <- carehomes_ifr_t_trajectories(steps, y, p,
                                       initial_step_from_parameters = FALSE)
  expect_equal(res1$step, matrix(steps, length(steps), ncol(y)))

  res2 <- carehomes_ifr_t_trajectories(steps, y, p,
                                       initial_step_from_parameters = TRUE)
  expect_equal(res2$step[1, ], step0)
  expect_equal(res2$step[-1, ], res1$step[-1, ])
})

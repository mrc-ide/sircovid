context("parameters")

test_that("single piecewise linear value", {
  expect_identical(sircovid_parameters_piecewise_linear(NULL, pi, 0.1), pi)
  expect_identical(sircovid_parameters_piecewise_linear(NULL,
                                                        matrix(c(1, pi),
                                                               nrow = 1),
                                                        0.1),
                   matrix(c(1, pi), nrow = 1))
  expect_error(sircovid_parameters_piecewise_linear(NULL, numeric(0), 0.1),
               "As 'date' is NULL, expected single value")
  expect_error(sircovid_parameters_piecewise_linear(NULL, 1:5, 0.1),
               "As 'date' is NULL, expected single value")
})


test_that("varying piecewise linear value", {
  ## TODO: this should get a better set of tests, as it's complex
  ## enough
  date <- as_sircovid_date(c("2020-02-01", "2020-02-10", "2020-02-29"))
  beta <- sircovid_parameters_piecewise_linear(date, 1:3, 0.5)
  expect_equal(
    beta,
    c(rep(1, 64),
      seq(1, 2, length.out = 19),
      seq(2, 3, length.out = 39)[-1]))

  date <- as_sircovid_date(c("2020-02-01", "2020-02-10", "2020-02-29"))
  beta <- sircovid_parameters_piecewise_linear(date, matrix(1:3, 3, 2), 0.5)
  expect_equal(
    beta,
    matrix(c(rep(1, 64),
      seq(1, 2, length.out = 19),
      seq(2, 3, length.out = 39)[-1]), 121, 2))
})


test_that("piecewise linear date and value have to be the same length", {
  date <- c(32, 41, 60)
  expect_error(
    sircovid_parameters_piecewise_linear(date, 1:2, 0.5),
    "'date' and 'value' must have the same length")
  expect_error(
    sircovid_parameters_piecewise_linear(date, 1:4, 0.5),
    "'date' and 'value' must have the same length")
})


test_that("can't use a single piecewise linear date/value", {
  expect_error(
    sircovid_parameters_piecewise_linear(32, 1, 0.5),
    "Need at least two dates and values for a varying piecewise linear")
})


test_that("piecewise linear dates must be increasing", {
  expect_error(
    sircovid_parameters_piecewise_linear(c(32, 41, 41, 64), 1:4, 0.5),
    "'date' must be strictly increasing")
})


test_that("piecewise linear dates must be sircovid_dates", {
  expect_error(
    sircovid_parameters_piecewise_linear(as_date(c("2020-02-01", "2020-02-10")),
                                         1:2, 0.5),
    "'date' must be numeric - did you forget sircovid_date()?")
  expect_error(
    sircovid_parameters_piecewise_linear(c(-10, 41, 60), 1:3, 0.5),
    "Negative dates, sircovid_date likely applied twice")
})


test_that("single piecewise constant value", {
  expect_identical(sircovid_parameters_piecewise_constant(NULL, pi, 0.1), pi)
  expect_error(sircovid_parameters_piecewise_constant(NULL, numeric(0), 0.1),
               "As 'date' is NULL, expected single value")
  expect_error(sircovid_parameters_piecewise_constant(NULL, 1:5, 0.1),
               "As 'date' is NULL, expected single value")
})


test_that("varying piecewise constant value", {
  date <- as_sircovid_date(c("2019-12-31", "2020-02-10", "2020-02-29"))
  y <- sircovid_parameters_piecewise_constant(date, 1:3, 0.5)
  expect_equal(
    y,
    c(rep(1, 82),
      rep(2, 38),
      3))
})


test_that("piecewise constant date and value have to be the same length", {
  date <- c(0, 41, 60)
  expect_error(
    sircovid_parameters_piecewise_constant(date, 1:2, 0.5),
    "'date' and 'value' must have the same length")
  expect_error(
    sircovid_parameters_piecewise_constant(date, 1:4, 0.5),
    "'date' and 'value' must have the same length")
})


test_that("piecewise constant first date must be 0", {
  expect_error(
    sircovid_parameters_piecewise_constant(c(20, 31, 41, 64), 1:4, 0.5),
    "As 'date' is not NULL, first date should be 0")
})


test_that("piecewise constant dates must be increasing", {
  expect_error(
    sircovid_parameters_piecewise_constant(c(0, 41, 41, 64), 1:4, 0.5),
    "'date' must be strictly increasing")
})


test_that("piecewise constant dates must be sircovid_dates", {
  expect_error(
    sircovid_parameters_piecewise_constant(
      as_date(c("2020-02-01", "2020-02-10")), 1:2, 0.5),
    "'date' must be numeric - did you forget sircovid_date()?")
  expect_error(
    sircovid_parameters_piecewise_constant(c(-10, 41, 60), 1:3, 0.5),
    "Negative dates, sircovid_date likely applied twice")
})


test_that("can read the default severity file", {
  data <- sircovid_parameters_severity(NULL)

  expect_identical(
    sircovid_parameters_severity(severity_default()),
    data)

  expect_vector_equal(lengths(data), 17)
  expect_setequal(
    names(data),
    c("p_star", "p_C", "p_G_D", "p_H_D",
      "p_ICU_D", "p_W_D", "p_ICU", "p_R",
      "p_sero_pos_1", "p_sero_pos_2", "p_H"))
  expect_vector_equal(data$p_serocoversion, data$p_serocoversion[[1]])
  expect_equal(
    data$p_G_D, rep(0.05, 17))
  expect_equal(
    data$p_star, rep(0.2, 17))
})


test_that("can validate a severity input", {
  d <- severity_default()
  expect_error(
    sircovid_parameters_severity(d[-1, ]),
    "Elements missing from 'data': 'p_C'")
})


test_that("can reprocess severity", {
  s <- sircovid_parameters_severity(NULL)
  expect_identical(
    sircovid_parameters_severity(s),
    s)
  expect_error(
    sircovid_parameters_severity(s[-1]),
    "Elements missing from 'params': 'p_star'")
})


test_that("shared parameters accepts a beta vector", {
  date <- sircovid_date("2020-02-01")
  beta_date <- sircovid_date(c("2020-02-01", "2020-02-14", "2020-03-15"))
  beta_value <- c(3, 1, 2)
  pars <- sircovid_parameters_shared(date, "england", beta_date, beta_value,
                                     "piecewise-linear", NULL, 1, 10)
  expect_equal(
    pars$beta_step,
    sircovid_parameters_piecewise_linear(beta_date, beta_value, 0.25))

  beta_date <- sircovid_date(c("2019-12-31", "2020-02-14", "2020-03-15"))
  beta_value <- c(3, 1, 2)
  pars <- sircovid_parameters_shared(date, "england", beta_date, beta_value,
                                     "piecewise-constant", NULL, 1, 10)
  expect_equal(
    pars$beta_step,
    sircovid_parameters_piecewise_constant(beta_date, beta_value, 0.25))

  expect_error(pars <- sircovid_parameters_shared(date, "england", beta_date,
                                                  beta_value,
                                                  "quadratic", NULL, 1, 10),
               "'beta_type' must be 'piecewise-linear' or 'piecewise-constant'")
})


test_that("shared parameters", {
  date <- sircovid_date("2020-02-01")
  pars <- sircovid_parameters_shared(date, "england", NULL, 0.1,
                                     "piecewise-linear", NULL, 1, 10)
  expect_setequal(
    names(pars),
    c("hosp_transmission", "ICU_transmission", "G_D_transmission",
      "dt", "steps_per_day", "n_age_groups",
      "beta_step", "population", "seed_step_start", "seed_value"))
  expect_equal(pars$beta_step, 0.1)
})


test_that("can expand beta", {
  date <- sircovid_date(c("2020-02-01", "2020-02-14", "2020-03-15"))
  value <- c(3, 1, 2)
  beta <- sircovid_parameters_piecewise_linear(date, value, 1)

  # The implied time series looks like this:
  t1 <- seq(0, date[[3]])
  res1 <- cbind(t1, beta, deparse.level = 0)

  expect_equal(sircovid_parameters_expand_step(t1, beta), beta)

  t2 <- seq(0, 100, by = 1)
  beta2 <- sircovid_parameters_expand_step(t2, beta)
  expect_equal(beta2[seq_along(beta)], beta)
  expect_equal(beta2[-seq_along(beta)], rep(beta[length(beta)], 25))

  t3 <- t2[1:65]
  beta3 <- sircovid_parameters_expand_step(t3, beta)
  expect_equal(beta3, beta[1:65])
})

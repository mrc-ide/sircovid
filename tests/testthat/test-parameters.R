context("parameters")

test_that("single beta value", {
  expect_identical(sircovid_parameters_beta(NULL, pi, 0.1), pi)
  expect_error(sircovid_parameters_beta(NULL, numeric(0), 0.1),
               "As 'date' is NULL, expected single value")
  expect_error(sircovid_parameters_beta(NULL, 1:5, 0.1),
               "As 'date' is NULL, expected single value")
})


test_that("varying beta value", {
  ## TODO: this should get a better set of tests, as it's complex
  ## enough
  date <- as_sircovid_date(c("2020-02-01", "2020-02-10", "2020-02-29"))
  beta <- sircovid_parameters_beta(date, 1:3, 0.5)
  expect_equal(
    beta,
    c(rep(1, 64),
      seq(1, 2, length.out = 19),
      seq(2, 3, length.out = 39)[-1]))
})


test_that("beta date and value have to be the same length", {
  date <- c(32, 41, 60)
  expect_error(
    sircovid_parameters_beta(date, 1:2, 0.5),
    "'date' and 'value' must have the same length")
  expect_error(
    sircovid_parameters_beta(date, 1:4, 0.5),
    "'date' and 'value' must have the same length")
})


test_that("can't use a single date/value", {
  expect_error(
    sircovid_parameters_beta(32, 1, 0.5),
    "Need at least two dates and values for a varying 'beta'")
})


test_that("dates must be increasing", {
  expect_error(
    sircovid_parameters_beta(c(32, 41, 41, 64), 1:4, 0.5),
    "'date' must be strictly increasing")
})


test_that("dates must be sircovid_dates", {
  expect_error(
    sircovid_parameters_beta(as_date(c("2020-02-01", "2020-02-10")), 1:2, 0.5),
    "'date' must be numeric - did you forget sircovid_date()?")
  expect_error(
    sircovid_parameters_beta(c(-10, 41, 60), 1:3, 0.5),
    "Negative dates, sircovid_date likely applied twice")
})


test_that("can read the default severity file", {
  data <- sircovid_parameters_severity(NULL)

  expect_identical(
    sircovid_parameters_severity(severity_default()),
    data)

  expect_true(all(lengths(data) == 17))
  expect_setequal(
    names(data),
    c("p_admit_conf", "p_asympt", "p_death_comm", "p_death_hosp_D",
      "p_death_ICU", "p_death_stepdown", "p_hosp_ILI", "p_ICU_hosp",
      "p_seroconversion", "p_sympt_ILI"))
  expect_true(
    all(data$p_serocoversion == data$p_serocoversion[[1]]))
  expect_equal(
    data$p_death_comm, rep(0, 17))
  expect_equal(
    data$p_admit_conf, rep(0.2, 17))
})


test_that("can validate a severity input", {
  d <- severity_default()
  expect_error(
    sircovid_parameters_severity(d[-2, ]),
    "Elements missing from 'data': 'Proportion with symptoms'")
})


test_that("can reprocess severity", {
  s <- sircovid_parameters_severity(NULL)
  expect_identical(
    sircovid_parameters_severity(s),
    s)
  expect_error(
    sircovid_parameters_severity(s[-1]),
    "Elements missing from 'params': 'p_admit_conf'")
})


test_that("shared parameters accepts a beta vector", {
  date <- sircovid_date("2020-02-01")
  beta_date <- sircovid_date(c("2020-02-01", "2020-02-14", "2020-03-15"))
  beta_value <- c(3, 1, 2)
  pars <- sircovid_parameters_shared(date, "england", beta_date, beta_value)
  expect_equal(
    pars$beta_step,
    sircovid_parameters_beta(beta_date, beta_value, 0.25))
})


test_that("shared parameters", {
  date <- sircovid_date("2020-02-01")
  pars <- sircovid_parameters_shared(date, "england", NULL, 0.1)
  expect_setequal(
    names(pars),
    c("hosp_transmission", "ICU_transmission", "comm_D_transmission",
      "dt", "initial_step", "n_age_groups", "beta_step", "population"))
  expect_equal(pars$beta_step, 0.1)
  expect_equal(pars$initial_step, date * 4)
})


test_that("can expand beta", {
  date <- sircovid_date(c("2020-02-01", "2020-02-14", "2020-03-15"))
  value <- c(3, 1, 2)
  beta <- sircovid_parameters_beta(date, value, 1)

  # The implied time series looks like this:
  t1 <- seq(0, date[[3]])
  res1 <- cbind(t1, beta, deparse.level = 0)

  expect_equal(sircovid_parameters_beta_expand(t1, beta), beta)

  t2 <- seq(0, 100, by = 1)
  beta2 <- sircovid_parameters_beta_expand(t2, beta)
  expect_equal(beta2[seq_along(beta)], beta)
  expect_equal(beta2[-seq_along(beta)], rep(beta[length(beta)], 25))

  t3 <- t2[1:65]
  beta3 <- sircovid_parameters_beta_expand(t3, beta)
  expect_equal(beta3, beta[1:65])
})

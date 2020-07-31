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
    "Need at least two dates and betas for a varying beta")
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
  expect_identical(
    sircovid_parameters_severity(severity_default()),
    data)

  expect_true(all(lengths(data) == 17))
  expect_setequal(
    names(data),
    c("p_asympt", "p_sympt_ILI", "p_recov_ICU", "p_recov_ILI", "p_hosp_ILI",
      "p_recov_hosp", "p_death_hosp", "p_death_hosp_D", "p_death_ICU",
      "p_ICU_hosp", "p_seroconversion", "p_death_comm", "p_admit_conf"))

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


test_that("shared parameters", {
  date <- sircovid_date("2020-02-01")
  pars <- sircovid_parameters_shared(date, "england", NULL, 0.1)
  expect_setequal(
    names(pars),
    c("hosp_transmission", "ICU_transmission", "comm_D_transmission",
      "dt", "initial_step", "N_age", "beta_step", "population"))
  expect_equal(pars$beta_step, 0.1)
})

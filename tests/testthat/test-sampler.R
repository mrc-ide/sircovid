context("sampler")

## This is not a real test but simply ties to ruin the model
test_that("sampler runs without error", {
  time_steps_per_day <- 4

  pars_model <- generate_parameter_rtm(
    seed_SSP=10,
    dt = 1/time_steps_per_day
  )

  pars_obs <- list(
    ## what should this be?
    phi_ICU = 0.95,
    ## what should this be?
    k_ICU = 2,
    ## current proportion of England deaths over UK deaths
    phi_death = 926/1019,
    ## what should this be?
    k_death = 2,
    #rate for exponential noise, something big so noise is small (but non-zero))
    exp_noise=1e6)

  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  data$date <- as.Date(data$date)

  pars_model$beta <- 0.125
  start_date <- as.Date("2020-02-02")
  X <- particle_filter(data = data, pars_model, pars_obs, n_particles = 100,
                       start_date = start_date,
                       time_steps_per_day = time_steps_per_day)

  expect_is(X, "list")
  expect_equal(names(X), "log_likelihood")
})

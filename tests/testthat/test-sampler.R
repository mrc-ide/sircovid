context("sampler")

## This is not a real test but simply tries to run the model
test_that("sampler runs without error", {
  set.seed(1)
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
  set.seed(1)
  X <- particle_filter(data = data, pars_model, pars_obs, n_particles = 100,
                       start_date = start_date,
                       time_steps_per_day = time_steps_per_day)

  expect_is(X, "list")
  expect_equal(names(X), "log_likelihood")

  set.seed(1)
  Y <- particle_filter(data = data, pars_model, pars_obs, n_particles = 100,
                       start_date = start_date,
                       time_steps_per_day = time_steps_per_day,
                       output_states = TRUE,
                       save_particles = TRUE)
  expect_equal(X$log_likelihood, Y$log_likelihood)
  expect_setequal(names(Y), c("log_likelihood", "states", "index"))
  ##                            t   state  particles
  expect_equal(dim(Y$states), c(58, 238,   100))

  ## saveRDS(Y, "reference.rds")
  expect_equal(Y, readRDS("reference.rds"))

  ## Testing plotting is always a nightmare
  if (FALSE) {
    date <- as.Date(start_date) + seq_len(nrow(Y$states)) - 1L
    particles <- apply(Y$states[, c(Y$index$I_ICU) - 1L, ], c(1, 3), sum)
    plot_particles(particles = particles,
                   particle_dates = date,
                   data = data$itu / pars_obs$phi_ICU,
                   data_dates = data$date,
                   ylab = "ICU")
  }
})

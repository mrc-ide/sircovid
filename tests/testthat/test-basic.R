context("basic")

test_that("can run the basic model", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england")
  mod <- basic$new(p, 0, 10)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- basic_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(basic_index(mod$info())$run)
  res <- mod$run(end)
  expected <-
    rbind(c(1352, 2102, 1714, 1498, 878, 729, 604, 1037, 1206, 1424),
          c(275842, 273785, 275310, 275167, 276831, 276762, 277848, 276574,
            276184, 275348))
  expect_equal(res, expected)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- basic_parameters(start_date, "england")
  data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                        start_date, pars$dt)

  n_particles <- 100
  pf <- mcstate::particle_filter$new(
    data, basic, n_particles,
    compare = basic_compare,
    initial = basic_initial,
    index = basic_index)

  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    # what should this be?
    phi_ICU = 0.95,
    # what should this be?
    k_ICU = 2,
    # current proportion of England deaths over UK deaths
    phi_death = 926 / 1019,
    # what should this be?
    k_death = 2,
    # rate for exponential noise, something big so noise is small (but
    # non-zero))
    exp_noise = 1e6)

  pf$run(pars, pars_obs, pars)
})

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
    rbind(c(1105, 1195, 1551, 2239, 5073, 1167, 1580, 1468, 1550, 2424),
          c(275845, 275199, 275025, 273361, 264499, 276074, 276170, 274926,
            275620, 272092))
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

  pars_obs <- basic_parameters_observation()

  ## TODO: this is not yet reliable, nor probably correct
  expect_type(pf$run(pars, pars_obs, pars), "double")
})

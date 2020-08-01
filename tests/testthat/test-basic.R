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
    rbind(c(1387, 2056, 1686, 1511, 933, 698, 551, 1072, 1140, 1501),
          c(274551, 274106, 275084, 276011, 277633, 277435, 278266, 275557,
            275412, 276108))
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

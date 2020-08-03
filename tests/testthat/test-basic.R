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
    rbind(c(81, 112, 313, 196, 228, 113, 101, 87, 60, 90),
          c(290040, 290250, 289270, 289768, 288659, 289774, 290519, 289045,
            290380, 289536))
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

  ## TODO: this is not yet reliable, nor probably correct
  ## TODO: mcstate should spread out parameters like this for us by default
  expect_type(pf$run(pars, pars, pars), "double")
})

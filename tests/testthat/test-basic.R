context("basic")

test_that("can run the basic model", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england")
  mod <- basic$new(p, 0, 10, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- basic_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(basic_index(mod$info())$run)
  res <- mod$run(end)

  expected <-
    rbind(icu = c(96, 23, 76, 67, 61, 81, 53, 59, 25, 103),
          deaths = c(267228, 268129, 269071, 267589, 267566, 267310, 268740,
                     267631, 267839, 267448))
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
  expect_type(pf$run(pars), "double")
})

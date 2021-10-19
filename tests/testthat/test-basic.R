context("basic")

test_that("can run the basic model", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england")
  mod <- basic$new(p, 0, 10, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- basic_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)

  mod$set_index(basic_index(mod$info())$run)
  res <- mod$run(end)

  ## Regnerate with: dput_named_matrix(res)
  expected <-
    rbind(icu        = c(0, 2, 1, 2, 0, 0, 1, 0, 0, 0),
          deaths     = c(310178, 310871, 310280, 308829, 309890, 309425,
                         310307, 310720, 310154, 309521),
          deaths_inc = c(0, 0, 1, 1, 1, 1, 0, 0, 0, 0))
  expect_equal(res, expected)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- basic_parameters(start_date, "england")
  data <- basic_data(read_csv(sircovid_file("extdata/example.csv")),
                     start_date, pars$dt)

  n_particles <- 100
  pf <- mcstate::particle_filter$new(
    data, basic, n_particles,
    compare = basic_compare,
    initial = basic_initial,
    index = basic_index)
  pf$run(pars)

  ## TODO: this is not yet reliable, nor probably correct
  ## TODO: mcstate should spread out parameters like this for us by default
  expect_type(pf$run(pars), "double")
})


test_that("incidence calculation is correct", {
  start_date <- sircovid_date("2020-02-02")
  pars <- basic_parameters(start_date, "england")
  mod <- basic$new(pars, 0, 10, n_threads = 10)
  info <- mod$info()
  initial <- basic_initial(info, 10, pars)
  mod$update_state(state = initial$state, step = initial$step)

  index <- c(D = info$index$D_tot,
             D_inc = info$index$D_inc)
  expect_length(index, 2) # guard against name changes

  steps <- seq(initial$step, length.out = 60 * 4 + 1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  i <- which(steps %% pars$steps_per_day == 0)
  y0 <- y[, , i[-length(i)]]
  y1 <- y[, , i[-1]]
  expect_equal(y1[1, , ] - y0[1, , ], y1[2, , ])
})


test_that("compiled compare function is correct", {
  start_date <- sircovid_date("2020-02-02")
  pars <- basic_parameters(start_date, "england", exp_noise = Inf)
  data <- basic_data(read_csv(sircovid_file("extdata/example.csv")),
                     start_date, pars$dt)

  np <- 10
  mod <- basic$new(pars, 0, np, seed = 1L)
  initial <- basic_initial(mod$info(), np, pars)
  mod$update_state(state = initial$state, step = initial$step)
  mod$set_index(basic_index(mod$info())$run)

  mod$set_data(dust::dust_data(data, "step_end"))

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  y <- mod$run(data$step_end[[i]])
  expect_equal(mod$compare_data(),
               basic_compare(y, data[i, ], pars))
})

context("basic")

test_that("can run the basic model", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england",
                        initial_seed_size = 10)
  mod <- basic$new(p, 0, 10, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- basic_initial(mod$info(), 10, p)
  mod$update_state(state = initial)

  mod$set_index(basic_index(mod$info())$run)
  res <- mod$run(end)

  ## Regnerate with: dput_named_matrix(res)
  expected <-
    rbind(icu               = c(1, 1, 0, 2, 0, 1, 1, 1, 1, 2),
          deaths            = c(309368, 309502, 311310, 309939, 310328,
                                309843, 310787, 309644, 310771, 308737),
          deaths_inc        = c(3, 2, 1, 0, 1, 0, 2, 0, 1, 1),
          fitted_icu        = c(0.950001049565452, 0.950000838137851,
                                3.13622893533116e-08, 1.90000099594009,
                                7.55899998800052e-07, 0.950000826396304,
                                0.950001256470541, 0.950000036892852,
                                0.950000094262866, 1.90000008517751),
          fitted_deaths_inc = c(2.72620268894599, 1.81746870631162,
                                0.908736431474323, 2.08020683292907e-06,
                                0.908735120406551, 9.43629704961186e-07,
                                1.81746895324511, 3.40637099819306e-07,
                                0.908734616454928, 0.908734236916313))
  expect_equal(res, expected)
})


test_that("initial seeding in one big lump", {
  start_date <- sircovid_date("2020-02-07")
  n_particles <- 20
  p <- basic_parameters(start_date, "england", initial_seed_size = 10,
                        initial_seed_pattern = NULL)
  mod <- basic$new(p, 4, n_particles, seed = 1L)
  end <- sircovid_date("2020-02-28") / p$dt

  initial <- basic_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]
  i <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  expect_equal(t[i],
               rep(start_date * p$steps_per_day + 1, n_particles))

  ## Total infections through seeding are plausible
  n <- mean(n_E[4, , i[[1]]])
  expect_gt(ppois(n, 10), 0.05)

  ## No natural infections in this period:
  expect_true(all(n_E[-4, , seq_len(i[[1]])] == 0))
})


test_that("initial seeding spread out", {
  start_date <- sircovid_date("2020-02-07") + 0.123
  n_particles <- 20
  pattern <- rep(1, 4) # over a 1 day window
  p <- basic_parameters(start_date, "england", initial_seed_size = 10,
                        initial_seed_pattern = pattern)

  expect_equal(p$seed_step_start, 152)
  expect_equal(p$seed_value, c(1.27, 2.5, 2.5, 2.5, 1.23))

  mod <- basic$new(p, 4, n_particles, seed = 1L)
  end <- sircovid_date("2020-02-28") / p$dt

  initial <- basic_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]
  i <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  expect_equal(min(t[i]), floor(start_date * p$steps_per_day) + 1)
  expect_gte(diff(range(t[i])), 1)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- basic_parameters(start_date, "england")
  data <- helper_basic_data(read_csv(sircovid_file("extdata/example.csv")),
                            0, pars$dt)

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
  mod$update_state(state = initial)

  index <- c(D = info$index$D_tot,
             D_inc = info$index$D_inc)
  expect_length(index, 2) # guard against name changes

  steps <- seq(0, length.out = 60 * 4 + 1)
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
  data <- helper_basic_data(read_csv(sircovid_file("extdata/example.csv")),
                            0, pars$dt)

  np <- 10
  mod <- basic$new(pars, 0, np, seed = 1L)
  initial <- basic_initial(mod$info(), np, pars)
  mod$update_state(state = initial)
  mod$set_index(basic_index(mod$info())$run)

  mod$set_data(dust::dust_data(data, "time_end"))

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  y <- mod$run(data$time_end[[i]])
  expect_equal(mod$compare_data(),
               basic_compare(y, data[i, ], pars))
})

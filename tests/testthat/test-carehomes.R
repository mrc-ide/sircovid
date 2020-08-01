context("carehomes")

test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)

  expected <- rbind(
    c(43,       8,        20,       39,       109),
    c(146,      52,       97,       112,      386),
    c(15844,    15862,    15751,    15942,    15561),
    c(295764,   296249,   295302,   294985,   296323),
    c(311608,   312111,   311051,   310926,   311878),
    c(145337,   145524,   146088,   145235,   145945),
    c(467527,   466021,   466788,   467135,   468575),
    c(6933,     2959,     4734,     5948,     16508),
    c(5245322,  5245255,  5243013,  5246940,  5240590),
    c(24086644, 24087335, 24089331, 24090745, 24077564))
  expect_equal(res, expected)
})


test_that("can run the particle filter on the model", {
  skip("Work in progress - needs new data")
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england")
  data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                        start_date, pars$dt)

  n_particles <- 100
  pf <- mcstate::particle_filter$new(
    data, carehomes, n_particles,
    compare = carehomes_compare,
    initial = carehomes_initial,
    index = carehomes_index)

  pars_obs <- carehomes_parameters_observation()

  pf$run(pars, pars_obs, pars)
})

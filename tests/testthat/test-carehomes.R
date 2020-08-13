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
    icu = c(4, 25, 2, 13, 5),
    general = c(31, 71, 11, 38, 34),
    deaths_comm = c(15825, 15971, 15955, 15914, 15882),
    deaths_hosp = c(305578, 305933, 305954, 305572, 306014),
    admitted = c(151349, 150757, 151137, 150569, 150964),
    new = c(483755, 484827, 484942, 483498, 483546),
    prob_pos = c(0.723911065300095, 0.723769622060999, 0.723852644013753,
                 0.723847861727815, 0.723785265474777))
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

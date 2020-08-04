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
    icu = c(4, 26, 2, 13, 6),
    general = c(32, 72, 11, 37, 33),
    deaths_comm = c(15825, 15971, 15955, 15914, 15882),
    deaths_hosp = c(305577, 305933, 305954, 305572, 306014),
    deaths_tot = c(321401, 321903, 321909, 321486, 321896),
    admitted = c(151349, 150757, 151137, 150569, 150964),
    new = c(483755, 484827, 484942, 483498, 483546),
    R_pre_15_64 = c(2033, 5320, 932, 3088, 2122),
    R_neg_15_64 = c(5470674, 5472544, 5474482, 5473012, 5475747),
    R_pos_15_64 = c(25149548, 25144489, 25147586, 25147430, 25145179))
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

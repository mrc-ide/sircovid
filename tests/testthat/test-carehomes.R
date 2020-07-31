context("carehomes")

## TODO: this test should come out and be replaced by something that
## does not load black-box parameters. This version just runs the model
test_that("can run the carehomes model", {
  d <- readRDS("carehomes.rds")
  mod <- carehomes$new(d$user, 0, 10)
  expect_equal(names(mod$index()), names(d$index)[-1])
  mod$set_state(d$state)
  mod$set_index(integer(0))
  res <- mod$run(300)
  expect_equal(res, matrix(numeric(0), 0, 10))
})


test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)

  expected <- rbind(
    c(37,       9,        34,       37,       108),
    c(160,      67,       95,       123,      432),
    c(15836,    15870,    15742,    15942,    15594),
    c(295120,   294381,   294951,   294656,   295499),
    c(310956,   310250,   310693,   310598,   311088),
    c(145515,   145252,   145714,   145480,   145361),
    c(467323,   466132,   467234,   466757,   466597),
    c(6639,     2793,     4689,     5971,     16401),
    c(5245130,  5244105,  5247884,  5241942,  5237992),
    c(24088623, 24093283, 24089393, 24085759, 24074383))
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

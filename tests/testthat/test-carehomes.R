context("carehomes")

test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)

  expected <- rbind(
    icu = c(8, 8, 3, 18, 6),
    general = c(52, 67, 20, 67, 19),
    deaths_comm = c(15882, 15928, 15628, 16120, 15893),
    deaths_hosp = c(305649, 305638, 305740, 306023, 306163),
    admitted = c(151360, 150612, 151197, 151101, 150386),
    new = c(484244, 483918, 483997, 484935, 484355),
    prob_pos = c(0.723860522215455, 0.723923883973212, 0.724058180618437,
                 0.723865988088634, 0.723862084297041))
  expect_equal(res, expected)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england")
  data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                        start_date, pars$dt)
  ## Add additional columns
  data$deaths_hosp <- data$deaths
  data$deaths_comm <- NA
  data$deaths <- NA
  data$general <- NA
  data$admitted <- NA
  data$new <- NA
  data$npos_15_64 <- NA
  data$ntot_15_64 <- NA

  pf <- carehomes_particle_filter(data, 10)
  expect_s3_class(pf, "particle_filter")

  pf$run(pars)
})

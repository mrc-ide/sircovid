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
    icu = c(2, 3, 21, 1, 2),
    general = c(18, 20, 84, 15, 18),
    deaths_comm = c(15889, 15844, 16062, 15937, 15812),
    deaths_hosp = c(311152, 311293, 310538, 312167, 311438),
    admitted = c(151368, 151094, 150770, 151905, 150889),
    new = c(484831, 484377, 484387, 484760, 484043),
    sero_pos = c(2747790, 2800948, 4110382, 2472164, 2717259),
    sympt_cases = c(30880920, 30877372, 30878183, 30879154, 30875434),
    sympt_cases_over25 = c(20922636, 20917672, 20918141, 20921444, 20916277),
    react_pos = c(382, 411, 1933, 248, 331))
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
  data$hosp <- NA
  data$admitted <- NA
  data$new <- NA
  data$new_admitted <- NA
  data$npos_15_64 <- NA
  data$ntot_15_64 <- NA
  data$pillar2_pos <- NA
  data$pillar2_tot <- NA
  data$pillar2_cases <- NA
  data$pillar2_over25_pos <- NA
  data$pillar2_over25_tot <- NA
  data$pillar2_over25_cases <- NA
  data$react_pos <- NA
  data$react_tot <- NA

  pf <- carehomes_particle_filter(data, 10)
  expect_s3_class(pf, "particle_filter")

  pf$run(pars)
})

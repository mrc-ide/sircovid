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
    icu = c(5, 6, 5, 6, 8),
    general = c(16, 25, 30, 72, 45),
    deaths_comm = c(15726, 15879, 16031, 15679, 15938),
    deaths_hosp = c(305424, 306065, 305690, 305730, 306131),
    admitted = c(150430, 151216, 150052, 151354, 151078),
    new = c(483598, 484219, 484135, 483636, 483624),
    sero_pos = c(2715747, 3158775, 2867125, 3827660, 3718908),
    sympt_cases = c(30873500, 30877735, 30866639, 30876245, 30874781),
    sympt_cases_over25 = c(20916377, 20921071, 20907224, 20918285, 20917518),
    react_pos = c(331, 623, 399, 1470, 1220))
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

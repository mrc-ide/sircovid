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
    icu = c(6, 3, 9, 8, 1),
    general = c(18, 13, 39, 29, 11),
    deaths_comm = c(23388, 23884, 23423, 23771, 23729),
    deaths_hosp = c(283558, 282782, 282692, 283285, 282627),
    admitted = c(131968, 131514, 131631, 132451, 130830),
    new = c(421562, 420698, 420256, 420628, 421371),
    sero_pos = c(2914674, 2728608, 3407428, 3363490, 2326639),
    sympt_cases = c(12978662, 12979632, 12985324, 12978410, 12969487),
    sympt_cases_over25 = c(10087808, 10088670, 10092882, 10087873, 10082966),
    react_pos = c(386, 301, 774, 632, 165))
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

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
    icu = c(7, 2, 4, 25, 1),
    general = c(45, 7, 31, 94, 4),
    deaths_comm = c(23709, 23366, 23496, 23473, 23558),
    deaths_hosp = c(282912, 283192, 283889, 282507, 282734),
    admitted = c(131816, 131736, 132251, 132228, 132547),
    new = c(419749, 420226, 421277, 420402, 420445),
    sero_pos = c(3505005, 2363427, 3403301, 4411116, 2079763),
    sympt_cases = c(12974263, 12972299, 12977695, 12974676, 12978585),
    sympt_cases_over25 = c(10081966, 10081915, 10089062, 10083360, 10086817),
    react_pos = c(840, 172, 739, 2129, 120))
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

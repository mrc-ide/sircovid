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
    icu = c(2, 6, 15, 6, 6),
    general = c(10, 27, 64, 37, 28),
    deaths_comm = c(23205, 23647, 23499, 23560, 23540),
    deaths_hosp = c(282901, 283505, 283869, 282925, 283674),
    admitted = c(132255, 132018, 131506, 132176, 132295),
    new = c(420402, 420940, 421167, 420670, 421032),
    sero_pos = c(2422335, 3293652, 4087610, 3435148, 3003784),
    sympt_cases = c(12979084, 12975265, 12981163, 12978228, 12979110),
    sympt_cases_over25 = c(10089623, 10084822, 10089091, 10089359, 10091136),
    react_pos = c(228, 688, 1562, 839, 443))
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

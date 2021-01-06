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
    icu = c(3, 4, 13, 0, 4),
    general = c(5, 29, 63, 7, 27),
    deaths_comm = c(23617, 23312, 23395, 23644, 23506),
    deaths_hosp = c(282682, 282661, 283958, 283492, 283718),
    admitted = c(131549, 132204, 132209, 131684, 131809),
    new = c(419996, 421195, 421543, 420917, 419875),
    sero_pos = c(2421282, 3290496, 3963011, 2281277, 3057878),
    sympt_cases = c(12974036, 12983365, 12978134, 12984141, 12974846),
    sympt_cases_over25 = c(10085961, 10094791, 10088378, 10093046, 10085132),
    react_pos = c(229, 638, 1488, 138, 499))
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

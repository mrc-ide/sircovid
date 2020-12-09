context("carehomes")

test_that("can run the carehomes model", {
  skip("reordered") # TODO: fix this after things settle
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)

  expected <- rbind(
    icu = c(9, 2, 8, 17, 0),
    general = c(32, 4, 28, 84, 2),
    deaths_comm = c(23422, 23610, 23647, 23328, 23435),
    deaths_hosp = c(283548, 282000, 283792, 282756, 283749),
    admitted = c(132264, 132345, 132206, 132051, 132407),
    new = c(421110, 419533, 420267, 421005, 420805),
    sero_pos = c(3505523, 2365680, 3398212, 4405999, 2082758),
    sympt_cases = c(12983933, 12981131, 12975658, 12974380, 12976576),
    sympt_cases_over25 = c(10091770, 10088347, 10086192, 10085069, 10086572),
    react_pos = c(817, 168, 709, 2104, 75))
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

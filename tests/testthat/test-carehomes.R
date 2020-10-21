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
    icu = c(7, 6, 10, 22, 1),
    general = c(12, 28, 21, 93, 9),
    deaths_comm = c(15755, 15490, 15810, 15957, 15742),
    deaths_hosp = c(306557, 305907, 305827, 306570, 306897),
    admitted = c(151071, 149943, 151206, 150846, 151831),
    new = c(484264, 485023, 484041, 484970, 484513),
    sero_pos = c(25153428, 25148482, 25147612, 25145074, 25152018),
    sympt_cases = c(30886202, 30881184, 30872615, 30878988, 30880678),
    sympt_cases_over25 = c(20922501, 20920145, 20913635, 20920082, 20921007),
    react_pos = c(311, 524, 415, 1906, 186))
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

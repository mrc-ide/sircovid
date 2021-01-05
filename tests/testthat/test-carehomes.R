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
    icu = c(5, 3, 7, 3, 1),
    general = c(30, 9, 36, 28, 5),
    deaths_comm = c(23519, 23695, 23618, 23507, 23382),
    deaths_hosp = c(282986, 283650, 282762, 282492, 283300),
    admitted = c(132159, 131898, 132171, 132141, 132265),
    new = c(420560, 420675, 419776, 420315, 420790),
    sero_pos = c(3274594, 2727065, 3288196, 3234381, 2321346),
    sympt_cases = c(12973798, 12978588, 12975921, 12975608, 12978116),
    sympt_cases_over25 = c(10085634, 10089781, 10085294, 10083958, 10087035),
    react_pos = c(673, 299, 691, 595, 171))
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

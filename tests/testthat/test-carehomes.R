context("carehomes")

test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)
  saveRDS(res,"res.rds")

  expected <- rbind(
    icu = c(6, 4, 25, 8, 3),
    general = c(15, 21, 86, 19, 9),
    deaths_comm = c(15899, 15792, 15652, 15646, 15794),
    deaths_hosp = c(306743, 305375, 305377, 305267, 306436),
    admitted = c(151004, 150346, 150990, 150638, 151293),
    new = c(485531, 484147, 483827, 484073, 483995),
    sero_pos = c(2749042, 2803987, 4084059, 2498915, 2720240),
    sympt_cases = c(30877047, 30868840, 30878810, 30874292, 30875140),
    sympt_cases_over25 = c(20921889, 20913052, 20919513, 20915348, 20915030),
    react_pos = c(387, 445, 1807, 309, 369))
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

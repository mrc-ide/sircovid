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
    icu = c(6, 6, 1, 41, 16),
    general = c(40, 32, 12, 161, 58),
    deaths_comm = c(15720, 15724, 15781, 15606, 15632),
    deaths_hosp = c(311136, 311394, 310985, 310533, 310945),
    admitted = c(150925, 151628, 151158, 150552, 150872),
    new = c(485636, 484836, 485070, 484621, 483675),
    sero_pos = c(3102098, 3273000, 2448111, 4788364, 3621341),
    sympt_cases = c(30879147, 30877432, 30875693, 30879699, 30876144),
    sympt_cases_over25 = c(20919326, 20919079, 20916216, 20918813, 20921402),
    react_pos = c(614, 738, 246, 3482, 1120))
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

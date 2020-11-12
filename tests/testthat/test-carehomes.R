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
    icu = c(4, 3, 11, 2, 6),
    general = c(22, 19, 88, 15, 16),
    deaths_comm = c(15854, 15951, 15857, 15922, 15953),
    deaths_hosp = c(310942, 310448, 311237, 311027, 311685),
    admitted = c(151233, 150302, 150173, 150836, 151722),
    new = c(483717, 484131, 484837, 485881, 484682),
    sero_pos = c(2743210, 2802715, 4066077, 2500065, 2716364),
    sympt_cases = c(30880046, 30880765, 30872087, 30885657, 30881724),
    sympt_cases_over25 = c(20920044, 20919886, 20918900, 20923143, 20922411),
    react_pos = c(405, 428, 1695, 299, 334))
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

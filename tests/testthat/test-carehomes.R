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
    icu = c(0, 8, 4, 18, 0),
    general = c(21, 42, 19, 78, 6),
    deaths_comm = c(23550, 23649, 23190, 23677, 23700),
    deaths_hosp = c(283153, 282636, 283277, 283282, 282841),
    admitted = c(132113, 131443, 132031, 132091, 132228),
    new = c(421160, 420378, 420748, 420486, 420306),
    sero_pos = c(2515967, 3623516, 3174151, 4214168, 2079285),
    sympt_cases = c(12979452, 12977017, 12981531, 12980423, 12981789),
    sympt_cases_over25 = c(10089990, 10087058, 10089856, 10091340, 10090948),
    react_pos = c(250, 1017, 616, 1735, 139))
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

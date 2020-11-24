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
    icu = c(2, 0, 0, 2, 1),
    general = c(8, 3, 8, 9, 2),
    deaths_comm = c(24175, 24056, 24158, 24142, 23854),
    deaths_hosp = c(548693, 548630, 549055, 548514, 548292),
    admitted = c(264464, 264696, 264692, 264548, 264479),
    new = c(848809, 848631, 848922, 848246, 848991),
    sero_pos = c(3103816, 2386911, 2859102, 3185104, 2094629),
    sympt_cases = c(34001062, 33998929, 34004500, 34000669, 33996176),
    sympt_cases_over25 = c(23368525, 23365378, 23372199, 23369888, 23363079),
    react_pos = c(98, 30, 74, 118, 19))
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

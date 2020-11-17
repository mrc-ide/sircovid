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
    icu = c(9, 1, 16, 31, 3),
    general = c(28, 7, 70, 125, 6),
    deaths_comm = c(15598, 16134, 15706, 15849, 15818),
    deaths_hosp = c(310973, 310989, 310119, 310553, 311027),
    admitted = c(150959, 150871, 150535, 150300, 151261),
    new = c(485079, 484652, 485073, 483866, 484580),
    sero_pos = c(3118709, 2091439, 3780081, 4555909, 2098412),
    sympt_cases = c(30873010, 30876219, 30875675, 30880916, 30876091),
    sympt_cases_over25 = c(20916388, 20917998, 20918774, 20919527, 20916743),
    react_pos = c(619, 145, 1243, 2836, 150))
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

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
    icu = c(6, 18, 4, 9, 0),
    general = c(40, 67, 16, 30, 11),
    deaths_comm = c(23216, 23741, 23509, 23427, 23406),
    deaths_hosp = c(283267, 283289, 283166, 282781, 283457),
    admitted = c(131463, 132034, 131945, 132491, 132296),
    new = c(420833, 420444, 420470, 420402, 421668),
    sero_pos = c(3641149, 4059301, 3005869, 3537502, 2313089),
    sympt_cases = c(12968858, 12970576, 12974640, 12975629, 12976166),
    sympt_cases_over25 = c(10081567, 10082673, 10083338, 10086722, 10089392),
    react_pos = c(924, 1497, 496, 857, 173))
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

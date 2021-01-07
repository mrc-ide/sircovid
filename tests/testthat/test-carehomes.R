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
    icu = c(1, 1, 11, 1, 5),
    general = c(15, 22, 67, 7, 17),
    deaths_comm = c(23605, 23640, 23209, 23201, 23777),
    deaths_hosp = c(283221, 282854, 283009, 283137, 283289),
    admitted = c(132445, 132085, 131936, 132476, 132280),
    new = c(421057, 419726, 420490, 420324, 421492),
    sero_pos = c(2418635, 3294346, 3968645, 2278172, 3056259),
    sympt_cases = c(12986034, 12976864, 12976260, 12977617, 12979726),
    sympt_cases_over25 = c(10094904, 10085897, 10086870, 10089267, 10087875),
    react_pos = c(204, 616, 1409, 169, 476))
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

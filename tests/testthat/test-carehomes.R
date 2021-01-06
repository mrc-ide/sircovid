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
    icu = c(1, 3, 13, 9, 2),
    general = c(18, 11, 32, 33, 6),
    deaths_comm = c(23461, 23455, 23511, 23602, 23831),
    deaths_hosp = c(283418, 282859, 283468, 283899, 283979),
    admitted = c(132062, 131092, 132417, 131941, 132222),
    new = c(420396, 420097, 420687, 420080, 420167),
    sero_pos = c(3047602, 2731490, 3400037, 3358751, 2331312),
    sympt_cases = c(12974166, 12979519, 12978951, 12969601, 12977552),
    sympt_cases_over25 = c(10085673, 10088913, 10087570, 10080936, 10086865),
    react_pos = c(561, 287, 771, 752, 165))
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

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
    icu = c(9, 8, 4, 5, 11),
    general = c(46, 19, 27, 33, 33),
    deaths_comm = c(23539, 23476, 23590, 23688, 23595),
    deaths_hosp = c(283024, 283197, 282842, 282603, 283856),
    admitted = c(131714, 132149, 131558, 132045, 131640),
    new = c(420569, 420959, 419979, 420255, 422025),
    sero_pos = c(3582598, 3075712, 2939204, 3394468, 3372980),
    sympt_cases = c(12973520, 12978403, 12974631, 12980995, 12978238),
    sympt_cases_over25 = c(10086076, 10088100, 10086706, 10091740, 10086705),
    react_pos = c(907, 502, 438, 710, 744))
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

context("carehomes")


test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  info <- mod$info()
  initial <- carehomes_initial(info, 10, p)
  mod$set_state(initial$state, initial$step)

  index <- c(carehomes_index(info)$run,
             deaths_carehomes = info$index[["D_carehomes_tot"]],
             deaths_comm = info$index[["D_comm_tot"]],
             deaths_hosp = info$index[["D_hosp_tot"]],
             admitted = info$index[["cum_admit_conf"]],
             diagnoses = info$index[["cum_new_conf"]],
             sympt_cases = info$index[["cum_sympt_cases"]],
             sympt_cases_over25 = info$index[["cum_sympt_cases_over25"]]
  )

  mod$set_index(index)
  res <- mod$run(end)

  ## Regenerate with: dput_named_matrix(res)
  expected <-
    rbind(time                               = c(213, 213, 213, 213,
                                                 213),
          icu                                = c(111, 299, 340, 59, 42),
          general                            = c(396, 1133, 1277, 303,
                                                 135),
          deaths_carehomes_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_inc                    = c(2, 6, 1, 1, 0),
          deaths_hosp_inc                    = c(29, 74, 85, 19, 6),
          admitted_inc                       = c(16, 10, 11, 5, 1),
          diagnoses_inc                      = c(24, 73, 76, 26, 5),
          sero_pos_1                         = c(5197368, 6910920, 7089785,
                                                 4716406, 3652078),
          sero_pos_2                         = c(5196882, 6913766, 7095050,
                                                 4717892, 3651559),
          sympt_cases_inc                    = c(294, 826, 918, 209,
                                                 86),
          sympt_cases_non_variant_inc        = c(294, 826, 918, 209,
                                                 86),
          sympt_cases_over25_inc             = c(247, 675, 757, 180,
                                                 65),
          sympt_cases_under15_inc            = c(28, 75, 87, 12, 11),
          sympt_cases_15_24_inc              = c(19, 76, 74, 17, 10),
          sympt_cases_25_49_inc              = c(95, 257, 267, 54, 21),
          sympt_cases_50_64_inc              = c(73, 204, 212, 59, 20),
          sympt_cases_65_79_inc              = c(47, 143, 164, 46, 15),
          sympt_cases_80_plus_inc            = c(32, 71, 114, 21, 9),
          sympt_cases_non_variant_over25_inc = c(247, 675, 757, 180,
                                                 65),
          react_pos                          = c(20534, 61844, 68919,
                                                 15064, 6328),
          deaths_carehomes                   = c(1683, 1663, 1639, 1705,
                                                 1722),
          deaths_comm                        = c(22849, 22627, 22434,
                                                 22845, 22594),
          deaths_hosp                        = c(198907, 197449, 198075,
                                                 198348, 197756),
          admitted                           = c(89732, 88887, 89875,
                                                 89461, 89478),
          diagnoses                          = c(297454, 296304, 297181,
                                                 297138, 297675),
          sympt_cases                        = c(9802437, 9799220, 9792545,
                                                 9799207, 9809255),
          sympt_cases_over25                 = c(7479954, 7475585, 7471847,
                                                 7475980, 7484963))
  expect_equal(res, expected)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england")
  data <- carehomes_data(read_csv(sircovid_file("extdata/example.csv")),
                         start_date, pars$dt)

  pf <- carehomes_particle_filter(data, 10)
  expect_s3_class(pf, "particle_filter")

  pf$run(pars)
})


test_that("incidence calculation is correct", {
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england")
  mod <- carehomes$new(pars, 0, 10, n_threads = 10)
  info <- mod$info()
  initial <- carehomes_initial(info, 10, pars)
  mod$set_state(initial$state, initial$step)

  ## We have interesting values by time 60, step 240
  ## There are 10 incidence variables, so we want to pull 20 variables
  index <- c(deaths_comm = info$index$D_comm_tot,
             deaths_comm_inc = info$index$D_comm_inc,
             deaths_carehomes = info$index$D_carehomes_tot,
             deaths_carehomes_inc = info$index$D_carehomes_inc,
             deaths_hosp = info$index$D_hosp_tot,
             deaths_hosp_inc = info$index$D_hosp_inc,
             admitted = info$index$cum_admit_conf,
             admitted_inc = info$index$admit_conf_inc,
             diagnoses = info$index$cum_new_conf,
             diagnoses_inc = info$index$new_conf_inc,
             sympt_cases = info$index$cum_sympt_cases,
             sympt_cases_inc = info$index$sympt_cases_inc,
             sympt_cases_over25 = info$index$cum_sympt_cases_over25,
             sympt_cases_over25_inc = info$index$sympt_cases_over25_inc,
             sympt_cases_non_variant =
               info$index$cum_sympt_cases_non_variant,
             sympt_cases_non_variant_inc =
               info$index$sympt_cases_non_variant_inc,
             sympt_cases_non_variant_over25 =
               info$index$cum_sympt_cases_non_variant_over25,
             sympt_cases_non_variant_over25_inc =
               info$index$sympt_cases_non_variant_over25_inc,
             infections = info$index$cum_infections,
             infections_inc = info$index$infections_inc)
  expect_length(index, 20) # guard against name changes

  steps <- seq(initial$step, length.out = 60 * 4 + 1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  i <- which(steps %% pars$steps_per_day == 0)
  j <- seq(1, 20, by = 2)
  y0 <- y[, , i[-length(i)]]
  y1 <- y[, , i[-1]]
  yd <- y1[j, , ] - y0[j, , ]
  yi <- y1[-j, , ]
  rownames(yd) <- rownames(yi) <- NULL
  expect_equal(yd, yi)
})


test_that("compiled compare function is correct", {
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england", exp_noise = Inf)
  data <- carehomes_data(read_csv(sircovid_file("extdata/example.csv")),
                         start_date, pars$dt)

  np <- 1
  mod <- carehomes$new(pars, 0, np, seed = 1L)
  initial <- carehomes_initial(mod$info(), np, pars)
  mod$set_state(initial$state, initial$step)
  mod$set_index(carehomes_index(mod$info())$run)

  mod$set_data(dust::dust_data(data, "step_end"))

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  y <- mod$run(data$step_end[[i]])
  expect_equal(mod$compare_data(),
               carehomes_compare(y, data[i, ], pars))
})


test_that("Test compiled carehomes components", {
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england", exp_noise = Inf)
  ## use a non-integer kappa for pillar 2 cases
  pars$kappa_pillar2_cases <- 2.5
  data <- carehomes_data(read_csv(sircovid_file("extdata/example.csv")),
                         start_date, pars$dt)

  np <- 10
  mod <- carehomes$new(pars, 0, np, seed = 1L)
  initial <- carehomes_initial(mod$info(), np, pars)
  mod$set_state(initial$state, initial$step)
  mod$set_index(carehomes_index(mod$info())$run)

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  step <- data$step_end[[i]]
  y <- mod$run(step)

  ## quickly bodge together some data:
  d <- data[i, ]
  d[setdiff(names(d), "step_end")] <- NA_real_
  update_data <- function(values) {
    d[names(values)] <- values
    d
  }

  partial <- list(
    c(icu = 50),
    c(general = 50),
    c(hosp = 50),
    c(deaths_hosp = 30, deaths_carehomes = 10, deaths_comm = 10),
    c(deaths_hosp = 30, deaths_non_hosp = 20),
    c(deaths = 50),
    c(admitted = 50),
    c(new_diagnoses = 50),
    c(new_all_admissions = 50),
    c(npos_15_64 = 10, ntot_15_64 = 50),
    c(pillar2_pos = 10, pillar2_tot = 50),
    c(pillar2_cases = 50),
    c(pillar2_over25_pos = 10, pillar2_over25_tot = 50),
    c(pillar2_over25_cases = 50),
    c(react_pos = 10, react_tot = 50),
    c(strain_non_variant = 10, strain_tot = 50))

  for (p in partial) {
    d_test <- update_data(p)
    mod$set_data(dust::dust_data(d_test, "step_end"))
    expect_equal(mod$compare_data(),
                 unname(carehomes_compare(y, d_test, pars)))
  }
})


test_that("can run the particle filter on the model 2", {
  skip_on_windows_gha()
  set.seed(2)
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england")
  data <- carehomes_data(read_csv(sircovid_file("extdata/example.csv")),
                         start_date, pars$dt)

  np <- 50
  pf1 <- carehomes_particle_filter(data, np, compiled_compare = FALSE,
                                   seed = 2)
  pf2 <- carehomes_particle_filter(data, np, compiled_compare = TRUE, seed = 2)
  ll1 <- pf1$run(pars)
  ll2 <- pf2$run(pars)
  expect_lt(abs(ll1 - ll2), 60)
})

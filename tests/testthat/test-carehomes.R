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
    rbind(icu                                = c(119, 160, 125, 711,
                                                 158),
          general                            = c(424, 690, 435, 2710,
                                                 621),
          deaths_carehomes_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_inc                    = c(2, 5, 0, 13, 1),
          deaths_hosp_inc                    = c(25, 26, 26, 156, 39),
          admitted_inc                       = c(10, 10, 5, 31, 10),
          diagnoses_inc                      = c(21, 42, 27, 144, 42),
          sero_pos                           = c(3623151, 4260384, 3725379,
                                                 6549934, 4178360),
          sympt_cases_inc                    = c(303, 463, 309, 1996,
                                                 440),
          sympt_cases_over25_inc             = c(240, 369, 255, 1638,
                                                 360),
          sympt_cases_non_variant_over25_inc = c(240, 369, 255, 1638,
                                                 360),
          react_pos                          = c(21612, 35419, 23463,
                                                 152554, 33518),
          deaths_carehomes                   = c(1717, 1646, 1738, 1731,
                                                 1698),
          deaths_comm                        = c(22426, 22605, 22226,
                                                 22322, 22522),
          deaths_hosp                        = c(196681, 196509, 196244,
                                                 194794, 196440),
          admitted                           = c(89188, 88901, 88840,
                                                 88355, 89090),
          diagnoses                          = c(295878, 294898, 295346,
                                                 293320, 295083),
          sympt_cases                        = c(9752024, 9748504, 9758086,
                                                 9749042, 9760996),
          sympt_cases_over25                 = c(7438927, 7435358, 7443739,
                                                 7435077, 7446244))
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
  ## There are 8 incidence variables, so we want to pull 16 variables
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
             sympt_cases_non_variant_over25 =
               info$index$cum_sympt_cases_non_variant_over25,
             sympt_cases_non_variant_over25_inc =
               info$index$sympt_cases_non_variant_over25_inc)
  expect_length(index, 16) # guard against name changes

  steps <- seq(initial$step, length.out = 60 * 4 + 1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  i <- which(steps %% pars$steps_per_day == 0)
  j <- seq(1, 16, by = 2)
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
  start_date <- sircovid_date("2020-02-02")
  pars <- carehomes_parameters(start_date, "england")
  data <- carehomes_data(read_csv(sircovid_file("extdata/example.csv")),
                         start_date, pars$dt)

  np <- 50
  pf1 <- carehomes_particle_filter(data, np, compiled_compare = FALSE)
  pf2 <- carehomes_particle_filter(data, np, compiled_compare = TRUE)
  ll1 <- pf1$run(pars)
  ll2 <- pf2$run(pars)
  expect_lt(abs(ll1 - ll2), 60)
})

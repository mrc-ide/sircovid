context("lancelot")


test_that("can run the lancelot model", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  mod <- lancelot$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  info <- mod$info()
  initial <- lancelot_initial(info, 10, p)
  mod$update_state(state = initial)

  index <- c(lancelot_index(info)$run,
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
          icu                                = c(102, 179, 49, 103, 61),
          general                            = c(366, 720, 187, 356,
                                                 246),
          deaths_carehomes_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_inc                    = c(2, 2, 0, 2, 3),
          deaths_hosp_inc                    = c(25, 34, 9, 16, 13),
          deaths_hosp_0_49_inc               = c(1, 0, 1, 1, 1),
          deaths_hosp_50_54_inc              = c(0, 1, 0, 1, 1),
          deaths_hosp_55_59_inc              = c(3, 3, 0, 0, 0),
          deaths_hosp_60_64_inc              = c(6, 4, 2, 2, 3),
          deaths_hosp_65_69_inc              = c(4, 12, 0, 3, 1),
          deaths_hosp_70_74_inc              = c(2, 4, 1, 2, 3),
          deaths_hosp_75_79_inc              = c(3, 1, 1, 1, 3),
          deaths_hosp_80_plus_inc            = c(6, 9, 4, 6, 1),
          admitted_inc                       = c(4, 8, 2, 4, 2),
          diagnoses_inc                      = c(21, 35, 13, 19, 12),
          sero_pos_1                         = c(5486055, 6511760, 4432346,
                                                 5349829, 5017412),
          sero_pos_2                         = c(5487698, 6510353, 4434378,
                                                 5349904, 5019052),
          sympt_cases_inc                    = c(244, 389, 114, 222,
                                                 140),
          sympt_cases_non_variant_inc        = c(244, 389, 114, 222,
                                                 140),
          sympt_cases_over25_inc             = c(205, 324, 103, 176,
                                                 119),
          sympt_cases_under15_inc            = c(24, 44, 5, 24, 12),
          sympt_cases_15_24_inc              = c(15, 21, 6, 22, 9),
          sympt_cases_25_49_inc              = c(65, 112, 38, 48, 39),
          sympt_cases_50_64_inc              = c(61, 99, 30, 56, 40),
          sympt_cases_65_79_inc              = c(41, 66, 19, 47, 24),
          sympt_cases_80_plus_inc            = c(38, 47, 16, 25, 16),
          sympt_cases_non_variant_over25_inc = c(205, 324, 103, 176,
                                                 119),
          react_pos                          = c(17991, 34665, 8306,
                                                 16312, 12742),
          deaths_carehomes                   = c(1578, 1710, 1631, 1667,
                                                 1638),
          deaths_comm                        = c(24416, 24268, 24286,
                                                 24320, 24418),
          deaths_hosp                        = c(211606, 211169, 211371,
                                                 211287, 212224),
          admitted                           = c(95681, 95928, 95499,
                                                 96117, 96310),
          diagnoses                          = c(316560, 316418, 316622,
                                                 316407, 315945),
          sympt_cases                        = c(10260026, 10257280,
                                                 10257985, 10262426, 10263935),
          sympt_cases_over25                 = c(7879490, 7882363, 7880146,
                                                 7882792, 7884763))
  expect_equal(res, expected)
})


test_that("initial seeding in one big lump", {
  start_date <- sircovid_date("2020-02-07")
  n_particles <- 20
  p <- lancelot_parameters(start_date, "england")
  mod <- lancelot$new(p, 4, n_particles, seed = 1L)
  end <- sircovid_date("2020-02-28") / p$dt

  initial <- lancelot_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]
  i <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  expect_equal(t[i],
               rep(start_date * p$steps_per_day + 1, n_particles))

  ## Total infections through seeding are plausible
  n <- mean(n_E[4, , i[[1]]])
  expect_gt(ppois(n, 10), 0.05)

  ## No natural infections in this period:
  expect_true(all(n_E[-4, , seq_len(i[[1]])] == 0))
})


test_that("initial seeding spread out", {
  start_date <- sircovid_date("2020-02-07") + 0.123
  n_particles <- 20
  pattern <- rep(1, 4) # over a 1 day window
  p <- lancelot_parameters(start_date, "england",
                           initial_seed_size = 10,
                           initial_seed_pattern = pattern)

  expect_equal(p$seed_step_start, 152)
  expect_equal(p$seed_value, c(1.27, 2.5, 2.5, 2.5, 1.23))

  mod <- lancelot$new(p, 4, n_particles, seed = 1L)
  end <- sircovid_date("2020-02-28") / p$dt

  initial <- lancelot_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]
  i <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  expect_equal(min(t[i]), floor(start_date * p$steps_per_day) + 1)
  expect_gte(diff(range(t[i])), 1)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england")
  data <- lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                        0, pars$dt)

  pf <- lancelot_particle_filter(data, 10)
  expect_s3_class(pf, "particle_filter")

  pf$run(pars)
})


test_that("incidence calculation is correct", {
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england")
  mod <- lancelot$new(pars, 0, 10, n_threads = 10)
  info <- mod$info()
  initial <- lancelot_initial(info, 10, pars)
  mod$update_state(state = initial)

  ## We have interesting values by time 60, step 240
  ## There are 10 incidence variables, so we want to pull 20 variables
  index <- c(deaths_comm = info$index$D_comm_tot,
             deaths_comm_inc = info$index$D_comm_inc,
             deaths_carehomes = info$index$D_carehomes_tot,
             deaths_carehomes_inc = info$index$D_carehomes_inc,
             deaths_hosp = info$index$D_hosp_tot,
             deaths_hosp_inc = info$index$D_hosp_inc,
             deaths_hosp_0_49 = info$index$D_hosp_0_49_tot,
             deaths_hosp_0_49_inc = info$index$D_hosp_0_49_inc,
             deaths_hosp_50_54 = info$index$D_hosp_50_54_tot,
             deaths_hosp_50_54_inc = info$index$D_hosp_50_54_inc,
             deaths_hosp_55_59 = info$index$D_hosp_55_59_tot,
             deaths_hosp_55_59_inc = info$index$D_hosp_55_59_inc,
             deaths_hosp_60_64 = info$index$D_hosp_60_64_tot,
             deaths_hosp_60_64_inc = info$index$D_hosp_60_64_inc,
             deaths_hosp_65_69 = info$index$D_hosp_65_69_tot,
             deaths_hosp_65_69_inc = info$index$D_hosp_65_69_inc,
             deaths_hosp_70_74 = info$index$D_hosp_70_74_tot,
             deaths_hosp_70_74_inc = info$index$D_hosp_70_74_inc,
             deaths_hosp_75_79 = info$index$D_hosp_75_79_tot,
             deaths_hosp_75_79_inc = info$index$D_hosp_75_79_inc,
             deaths_hosp_80_plus = info$index$D_hosp_80_plus_tot,
             deaths_hosp_80_plus_inc = info$index$D_hosp_80_plus_inc,
             deaths = info$index$D_tot,
             deaths_inc = info$index$D_inc,
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
  expect_length(index, 38) # guard against name changes

  steps <- seq(0, length.out = 60 * 4 + 1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  i <- which(steps %% pars$steps_per_day == 0)
  j <- seq(1, length(index), by = 2)
  y0 <- y[, , i[-length(i)]]
  y1 <- y[, , i[-1]]
  yd <- y1[j, , ] - y0[j, , ]
  yi <- y1[-j, , ]
  rownames(yd) <- rownames(yi) <- NULL
  expect_equal(yd, yi)
})


test_that("compiled compare function is correct", {
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england", exp_noise = Inf)
  data <- lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                        0, pars$dt)

  np <- 1
  mod <- lancelot$new(pars, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, pars)
  mod$update_state(state = initial)
  mod$set_index(lancelot_index(mod$info())$run)

  mod$set_data(dust::dust_data(data, "step_end"))

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  y <- mod$run(data$step_end[[i]])
  expect_equal(mod$compare_data(),
               lancelot_compare(y, data[i, ], pars))
})


test_that("Test compiled lancelot components", {
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england", exp_noise = Inf)
  ## use a non-integer kappa for pillar 2 cases
  pars$kappa_pillar2_cases <- 2.5
  data <- lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                        0, pars$dt)

  np <- 10
  mod <- lancelot$new(pars, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, pars)
  mod$update_state(state = initial)
  mod$set_index(lancelot_index(mod$info())$run)

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
    c(deaths_hosp_0_49 = 1),
    c(deaths_hosp_50_54 = 1),
    c(deaths_hosp_55_59 = 1),
    c(deaths_hosp_60_64 = 2),
    c(deaths_hosp_65_69 = 2),
    c(deaths_hosp_70_74 = 3),
    c(deaths_hosp_75_79 = 5),
    c(deaths_hosp_80_plus = 10),
    c(deaths_hosp = 30),
    c(deaths_carehomes = 10),
    c(deaths_comm = 10),
    c(deaths_non_hosp = 20),
    c(deaths_hosp = 30, deaths_carehomes = 10, deaths_comm = 10),
    c(deaths_hosp = 30, deaths_non_hosp = 20),
    c(deaths = 50),
    c(admitted = 50),
    c(new_diagnoses = 50),
    c(new_all_admissions = 50),
    c(npos_15_64 = 10, ntot_15_64 = 50),
    c(pillar2_pos = 10, pillar2_tot = 50),
    c(pillar2_over25_pos = 10, pillar2_over25_tot = 50),
    c(pillar2_under15_pos = 10, pillar2_under15_tot = 50),
    c(pillar2_15_24_pos = 10, pillar2_15_24_tot = 50),
    c(pillar2_25_49_pos = 10, pillar2_25_49_tot = 50),
    c(pillar2_50_64_pos = 10, pillar2_50_64_tot = 50),
    c(pillar2_65_79_pos = 10, pillar2_65_79_tot = 50),
    c(pillar2_80_plus_pos = 10, pillar2_80_plus_tot = 50),
    c(pillar2_cases = 50),
    c(pillar2_over25_cases = 50),
    c(pillar2_under15_cases = 50),
    c(pillar2_15_24_cases = 50),
    c(pillar2_25_49_cases = 50),
    c(pillar2_50_64_cases = 50),
    c(pillar2_65_79_cases = 50),
    c(pillar2_80_plus_cases = 50),
    c(react_pos = 10, react_tot = 50),
    c(strain_non_variant = 10, strain_tot = 50))

  for (p in partial) {
    d_test <- update_data(p)
    mod$set_data(dust::dust_data(d_test, "step_end"))
    expect_equal(mod$compare_data(),
                 unname(lancelot_compare(y, d_test, pars)))
  }
})


test_that("can run the particle filter on the model 2", {
  skip_on_windows_gha()
  set.seed(2)
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england")
  data <- lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                        0, pars$dt)

  np <- 50
  pf1 <- lancelot_particle_filter(data, np, compiled_compare = FALSE,
                                  seed = 2)
  pf2 <- lancelot_particle_filter(data, np, compiled_compare = TRUE, seed = 2)
  ll1 <- pf1$run(pars)
  ll2 <- pf2$run(pars)
  expect_lt(abs(ll1 - ll2), 60)
})

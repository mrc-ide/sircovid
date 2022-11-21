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
          icu                                = c(23, 94, 62, 41, 72),
          general                            = c(141, 344, 320, 184,
                                                 325),
          deaths_carehomes_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_inc                    = c(1, 3, 2, 0, 3),
          deaths_comm_0_49_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_50_54_inc              = c(0, 0, 0, 0, 0),
          deaths_comm_55_59_inc              = c(0, 0, 0, 0, 0),
          deaths_comm_60_64_inc              = c(0, 1, 1, 0, 0),
          deaths_comm_65_69_inc              = c(0, 0, 0, 0, 2),
          deaths_comm_70_74_inc              = c(0, 1, 0, 0, 0),
          deaths_comm_75_79_inc              = c(0, 0, 0, 0, 0),
          deaths_comm_80_plus_inc            = c(1, 1, 1, 0, 1),
          deaths_hosp_inc                    = c(8, 21, 9, 10, 20),
          deaths_hosp_0_49_inc               = c(0, 2, 0, 0, 0),
          deaths_hosp_50_54_inc              = c(0, 1, 0, 1, 1),
          deaths_hosp_55_59_inc              = c(2, 0, 0, 0, 4),
          deaths_hosp_60_64_inc              = c(2, 4, 1, 3, 2),
          deaths_hosp_65_69_inc              = c(0, 2, 2, 1, 2),
          deaths_hosp_70_74_inc              = c(1, 4, 1, 2, 6),
          deaths_hosp_75_79_inc              = c(0, 4, 3, 1, 1),
          deaths_hosp_80_plus_inc            = c(3, 4, 2, 2, 4),
          admitted_inc                       = c(1, 2, 2, 2, 4),
          all_admission_0_9_inc              = c(0, 0, 0, 0, 0),
          all_admission_10_19_inc            = c(0, 0, 0, 0, 0),
          all_admission_20_29_inc            = c(0, 0, 0, 0, 0),
          all_admission_30_39_inc            = c(0, 0, 0, 1, 1),
          all_admission_40_49_inc            = c(0, 1, 1, 0, 1),
          all_admission_50_59_inc            = c(3, 1, 2, 0, 6),
          all_admission_60_69_inc            = c(1, 4, 3, 3, 5),
          all_admission_70_79_inc            = c(2, 2, 4, 3, 4),
          all_admission_80_plus_inc          = c(2, 7, 10, 0, 4),
          diagnoses_inc                      = c(7, 13, 18, 5, 17),
          sero_pos_1                         = c(4070155, 5315186, 5258910,
                                                 4437996, 5151040),
          sero_pos_2                         = c(4072810, 5313592, 5260266,
                                                 4435658, 5149290),
          sympt_cases_inc                    = c(67, 213, 188, 98, 181),
          sympt_cases_non_variant_inc        = c(67, 213, 188, 98, 181),
          sympt_cases_over25_inc             = c(51, 169, 163, 87, 150),
          sympt_cases_under15_inc            = c(7, 31, 9, 7, 19),
          sympt_cases_15_24_inc              = c(9, 13, 16, 4, 12),
          sympt_cases_25_49_inc              = c(15, 62, 59, 28, 59),
          sympt_cases_50_64_inc              = c(19, 51, 48, 26, 48),
          sympt_cases_65_79_inc              = c(10, 34, 35, 20, 31),
          sympt_cases_80_plus_inc            = c(7, 22, 21, 13, 12),
          sympt_cases_non_variant_over25_inc = c(51, 169, 163, 87, 150),
          ons_pos                            = c(6164.2, 16174.4, 15933.4,
                                                 8580.8, 14639.6),
          react_pos                          = c(5938, 15575, 15361,
                                                 8255, 14093),
          react_5_24_pos                     = c(1316, 3254, 3198, 1751,
                                                 2920),
          react_25_34_pos                    = c(851, 2215, 2291, 1229.25,
                                                 2024.25),
          react_35_44_pos                    = c(758, 2006, 2030, 1084.25,
                                                 1852.25),
          react_45_54_pos                    = c(828, 2393, 2353, 1244.25,
                                                 2127.25),
          react_55_64_pos                    = c(922, 2394, 2294, 1232.25,
                                                 2106.25),
          react_65_plus_pos                  = c(1263, 3313, 3195, 1714,
                                                 3063),
          deaths_carehomes                   = c(1635, 1624, 1639, 1727,
                                                 1659),
          deaths_comm                        = c(24455, 24379, 24237,
                                                 23976, 24328),
          deaths_hosp                        = c(211893, 211891, 211440,
                                                 211329, 211850),
          admitted                           = c(96462, 96241, 95901,
                                                 96567, 96072),
          diagnoses                          = c(316693, 316674, 316747,
                                                 316754, 316549),
          sympt_cases                        = c(10269626, 10270162,
                                                 10261809, 10256631, 10268512),
          sympt_cases_over25                 = c(7889147, 7890565, 7881589,
                                                 7877369, 7888184))
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
  data <- helper_lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                               0, pars$dt)

  pf <- helper_lancelot_particle_filter(data, 10)
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
  data <- helper_lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                               0, pars$dt)

  np <- 1
  mod <- lancelot$new(pars, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, pars)
  mod$update_state(state = initial)
  mod$set_index(lancelot_index(mod$info())$run)

  mod$set_data(dust::dust_data(data, "time_end"))

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  y <- mod$run(data$time_end[[i]])
  expect_equal(mod$compare_data(),
               lancelot_compare(y, data[i, ], pars))
})


test_that("Test compiled lancelot components", {
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england", exp_noise = Inf)
  ## use a non-integer kappa for pillar 2 cases
  pars$kappa_pillar2_cases <- 2.5
  data <- helper_lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                               0, pars$dt)

  np <- 10
  mod <- lancelot$new(pars, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, pars)
  mod$update_state(state = initial)
  mod$set_index(lancelot_index(mod$info())$run)

  i <- which(!is.na(data$icu) & !is.na(data$deaths))[[10]]
  time <- data$time_end[[i]]
  y <- mod$run(time)

  ## quickly bodge together some data:
  d <- data[i, ]
  d[setdiff(names(d), "time_end")] <- NA_real_
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
    c(deaths_comm_0_49 = 1),
    c(deaths_comm_50_54 = 1),
    c(deaths_comm_55_59 = 1),
    c(deaths_comm_60_64 = 2),
    c(deaths_comm_65_69 = 2),
    c(deaths_comm_70_74 = 3),
    c(deaths_comm_75_79 = 5),
    c(deaths_comm_80_plus = 10),
    c(deaths_comm = 10),
    c(deaths_non_hosp = 20),
    c(deaths = 50),
    c(admitted = 50),
    c(new_diagnoses = 50),
    c(new_all_admissions = 50),
    c(all_admission_0_9_inc = 0),
    c(all_admission_10_19_inc = 0),
    c(all_admission_20_29_inc = 2),
    c(all_admission_30_39_inc = 1),
    c(all_admission_40_49_inc = 0),
    c(all_admission_50_59_inc = 4),
    c(all_admission_60_69_inc = 5),
    c(all_admission_70_79_inc = 5),
    c(all_admission_80_plus_inc = 8),
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
    c(ons_pos = 20, ons_tot = 60),
    c(react_pos = 10, react_tot = 50),
    c(react_5_24_pos = 1, react_5_24_tot = 5),
    c(react_25_34_pos = 1, react_25_34_tot = 5),
    c(react_35_44_pos = 2, react_35_44_tot = 10),
    c(react_45_54_pos = 2, react_45_54_tot = 10),
    c(react_55_64_pos = 2, react_55_64_tot = 10),
    c(react_65_plus_pos = 2, react_65_plus_tot = 10),
    c(strain_non_variant = 10, strain_tot = 50))

  for (p in partial) {
    d_test <- update_data(p)
    mod$set_data(dust::dust_data(d_test, "time_end"))
    expect_equal(mod$compare_data(),
                 unname(lancelot_compare(y, d_test, pars)))
  }
})


test_that("can run the particle filter on the model 2", {
  skip_on_windows_gha()
  set.seed(2)
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england")
  data <- helper_lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                               0, pars$dt)

  np <- 50
  pf1 <- helper_lancelot_particle_filter(data, np, compiled_compare = FALSE,
                                         seed = 2)
  pf2 <- helper_lancelot_particle_filter(data, np, compiled_compare = TRUE,
                                         seed = 2)
  ll1 <- pf1$run(pars)
  ll2 <- pf2$run(pars)
  expect_lt(abs(ll1 - ll2), 60)
})

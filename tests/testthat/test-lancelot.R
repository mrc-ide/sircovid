context("lancelot")


test_that("can run the lancelot model", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  mod <- lancelot$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  info <- mod$info()
  initial <- lancelot_initial(info, 10, p)
  mod$update_state(state = initial$state, step = initial$step)

  index <- c(lancelot_index(info)$run,
             deaths_carehomes = info$index[["D_carehomes_tot"]],
             deaths_comm = info$index[["D_comm_tot"]],
             deaths_hosp = info$index[["D_hosp_tot"]],
             all_admission = info$index[["cum_all_admission"]],
             sympt_cases = info$index[["cum_sympt_cases"]],
             sympt_cases_over25 = info$index[["cum_sympt_cases_over25"]]
  )

  mod$set_index(index)
  res <- mod$run(end)

  ## Regenerate with: dput_named_matrix(res)
  expected <-
    rbind(time                               = c(213, 213, 213, 213,
                                                 213),
          icu                                = c(283, 333, 260, 290,
                                                 185),
          general                            = c(1135, 1269, 1130, 1042,
                                                 720),
          deaths_carehomes_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_inc                    = c(4, 6, 3, 7, 6),
          deaths_hosp_inc                    = c(52, 61, 45, 46, 35),
          all_admission_inc                  = c(64, 63, 51, 57, 39),
          sero_pos_1                         = c(7097636, 7402837, 7073092,
                                                 7037479, 6454816),
          sero_pos_2                         = c(7096365, 7400681, 7073068,
                                                 7032728, 6455180),
          sympt_cases_inc                    = c(634, 706, 618, 573,
                                                 430),
          sympt_cases_non_variant_inc        = c(634, 706, 618, 573,
                                                 430),
          sympt_cases_over25_inc             = c(521, 586, 535, 483,
                                                 355),
          sympt_cases_under15_inc            = c(65, 71, 50, 56, 41),
          sympt_cases_15_24_inc              = c(48, 49, 33, 34, 34),
          sympt_cases_25_49_inc              = c(190, 209, 195, 157,
                                                 111),
          sympt_cases_50_64_inc              = c(152, 172, 144, 153,
                                                 120),
          sympt_cases_65_79_inc              = c(109, 130, 133, 117,
                                                 78),
          sympt_cases_80_plus_inc            = c(70, 75, 63, 56, 46),
          sympt_cases_non_variant_over25_inc = c(521, 586, 535, 483,
                                                 355),
          react_pos                          = c(49406, 57931, 47989,
                                                 46491, 32864),
          deaths_carehomes                   = c(1754, 1687, 1707, 1637,
                                                 1663),
          deaths_comm                        = c(24190, 24491, 24586,
                                                 24281, 24310),
          deaths_hosp                        = c(211169, 210808, 211365,
                                                 212195, 211612),
          all_admission                      = c(492974, 491650, 493114,
                                                 494621, 493315),
          sympt_cases                        = c(10254718, 10244341,
                                                 10255725, 10268529, 10264053),
          sympt_cases_over25                 = c(7879262, 7869797, 7878007,
                                                 7889211, 7883259))

  expect_equal(res, expected)
})


test_that("can run the particle filter on the model", {
  start_date <- sircovid_date("2020-02-02")
  pars <- lancelot_parameters(start_date, "england")
  data <- lancelot_data(read_csv(sircovid_file("extdata/example.csv")),
                        start_date, pars$dt)

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
  mod$update_state(state = initial$state, step = initial$step)

  ## We have interesting values by time 60, step 240
  ## There are 10 incidence variables, so we want to pull 20 variables
  index <- c(deaths_comm = info$index$D_comm_tot,
             deaths_comm_inc = info$index$D_comm_inc,
             deaths_carehomes = info$index$D_carehomes_tot,
             deaths_carehomes_inc = info$index$D_carehomes_inc,
             deaths_hosp = info$index$D_hosp_tot,
             deaths_hosp_inc = info$index$D_hosp_inc,
             all_admission = info$index$cum_all_admission,
             all_admission_inc = info$index$all_admission_inc,
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
  expect_length(index, 18) # guard against name changes

  steps <- seq(initial$step, length.out = 60 * 4 + 1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  i <- which(steps %% pars$steps_per_day == 0)
  j <- seq(1, 18, by = 2)
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
                        start_date, pars$dt)

  np <- 1
  mod <- lancelot$new(pars, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, pars)
  mod$update_state(state = initial$state, step = initial$step)
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
                        start_date, pars$dt)

  np <- 10
  mod <- lancelot$new(pars, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, pars)
  mod$update_state(state = initial$state, step = initial$step)
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
    c(deaths_hosp = 30, deaths_carehomes = 10, deaths_comm = 10),
    c(deaths_hosp = 30, deaths_non_hosp = 20),
    c(deaths = 50),
    c(all_admission = 50),
    c(sero_pos_15_64_1 = 10, sero_tot_15_64_1 = 50),
    c(sero_pos_15_64_2 = 10, sero_tot_15_64_2 = 50),
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
                        start_date, pars$dt)

  np <- 50
  pf1 <- lancelot_particle_filter(data, np, compiled_compare = FALSE,
                                  seed = 2)
  pf2 <- lancelot_particle_filter(data, np, compiled_compare = TRUE, seed = 2)
  ll1 <- pf1$run(pars)
  ll2 <- pf2$run(pars)
  expect_lt(abs(ll1 - ll2), 60)
})

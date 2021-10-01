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
          icu                                = c(114, 400, 272, 76, 45),
          general                            = c(442, 1756, 985, 337,
                                                 134),
          deaths_carehomes_inc               = c(0, 0, 0, 0, 0),
          deaths_comm_inc                    = c(3, 9, 0, 0, 0),
          deaths_hosp_inc                    = c(31, 77, 49, 23, 8),
          all_admission_inc                  = c(34, 123, 63, 22, 12),
          sero_pos_1                         = c(5196621, 7338502, 6426135,
                                                 4706995, 3632879),
          sero_pos_2                         = c(5197954, 7336944, 6432145,
                                                 4704353, 3633934),
          sympt_cases_inc                    = c(266, 992, 571, 181,
                                                 77),
          sympt_cases_non_variant_inc        = c(266, 992, 571, 181,
                                                 77),
          sympt_cases_over25_inc             = c(215, 805, 476, 153,
                                                 65),
          sympt_cases_under15_inc            = c(28, 113, 55, 15, 9),
          sympt_cases_15_24_inc              = c(23, 74, 40, 13, 3),
          sympt_cases_25_49_inc              = c(78, 271, 164, 62, 17),
          sympt_cases_50_64_inc              = c(60, 244, 135, 45, 23),
          sympt_cases_65_79_inc              = c(51, 165, 112, 30, 15),
          sympt_cases_80_plus_inc            = c(26, 125, 65, 16, 10),
          sympt_cases_non_variant_over25_inc = c(215, 805, 476, 153,
                                                 65),
          react_pos                          = c(21039, 79104, 45801,
                                                 14807, 6286),
          deaths_carehomes                   = c(1735, 1754, 1609, 1617,
                                                 1706),
          deaths_comm                        = c(22849, 22500, 22680,
                                                 22528, 22787),
          deaths_hosp                        = c(197777, 198305, 198223,
                                                 198552, 198178),
          all_admission                      = c(463627, 463990, 463660,
                                                 464490, 463641),
          sympt_cases                        = c(9807206, 9802506, 9794172,
                                                 9803463, 9813121),
          sympt_cases_over25                 = c(7483748, 7480192, 7472728,
                                                 7481975, 7485904))
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

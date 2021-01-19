context("carehomes")

test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  info <- mod$info()
  initial <- carehomes_initial(info, 10, p)
  mod$set_state(initial$state, initial$step)
  
  index <- c(carehomes_index(info)$run,
             deaths_comm = info$index[["D_comm_tot"]],
             deaths_hosp = info$index[["D_hosp_tot"]],
             admitted = info$index[["cum_admit_conf"]],
             new = info$index[["cum_new_conf"]],
             sympt_cases = info$index[["cum_sympt_cases"]],
             sympt_cases_over25 = info$index[["cum_sympt_cases_over25"]]
  )
  
  mod$set_index(index)
  res <- mod$run(end)

  ## Regenerate with: dput_named_matrix(res)
  expected <-
    rbind(icu                                = c(3, 4, 5, 10, 3),
          general                            = c(35, 34, 35, 41, 7),
          deaths_comm_inc                    = c(0, 0, 0, 0, 0),
          deaths_hosp_inc                    = c(2, 1, 3, 2, 0),
          admitted_inc                       = c(0, 0, 2, 0, 0),
          new_inc                            = c(1, 1, 1, 2, 0),
          sero_pos                           = c(3043375, 3303884, 3660340,
                                                 3555769, 2384929),
          sympt_cases_inc                    = c(6, 3, 5, 8, 5),
          sympt_cases_over25_inc             = c(5, 3, 4, 6, 5),
          sympt_cases_non_variant_over25_inc = c(5, 3, 4, 6, 5),
          react_pos                          = c(470, 641, 1063, 877,
                                                 151),
          deaths_comm                        = c(23349, 23597, 23357,
                                                 23379, 23245),
          deaths_hosp                        = c(283383, 283433, 282859,
                                                 282832, 284127),
          admitted                           = c(131981, 132378, 132257,
                                                 131758, 132269),
          new                                = c(420336, 421268, 419989,
                                                 420711, 421502),
          sympt_cases                        = c(12977835, 12976953,
                                                 12983125, 12974885, 12981657),
          sympt_cases_over25                 = c(10087409, 10086355,
                                                 10092710, 10087814, 10091423))
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
  data$strain_non_variant <- NA
  data$strain_tot <- NA

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
  ## There are 7 incidence variables, so we want to pull 14 variables
  index <- c(deaths_comm = info$index$D_comm_tot,
             deaths_comm_inc = info$index$D_comm_inc,
             deaths_hosp = info$index$D_hosp_tot,
             deaths_hosp_inc = info$index$D_hosp_inc,
             admitted = info$index$cum_admit_conf,
             admitted_inc = info$index$admit_conf_inc,
             new = info$index$cum_new_conf,
             new_inc = info$index$new_conf_inc,
             sympt_cases = info$index$cum_sympt_cases,
             sympt_cases_inc = info$index$sympt_cases_inc,
             sympt_cases_over25 = info$index$cum_sympt_cases_over25,
             sympt_cases_over25_inc = info$index$sympt_cases_over25_inc,
             sympt_cases_non_variant_over25 =
               info$index$cum_sympt_cases_non_variant_over25,
             sympt_cases_non_variant_over25_inc =
               info$index$sympt_cases_non_variant_over25_inc)
  expect_length(index, 14) # guard against name changes

  steps <- seq(initial$step, length.out = 60 * 4 + 1)
  y <- dust::dust_iterate(mod, steps, index)

  i <- which(steps %% pars$steps_per_day == 0)
  j <- seq(1, 14, by = 2)
  y0 <- y[, , i[-length(i)]]
  y1 <- y[, , i[-1]]
  yd <- y1[j, , ] - y0[j, , ]
  yi <- y1[-j, , ]
  rownames(yd) <- rownames(yi) <- NULL
  expect_equal(yd, yi)
})

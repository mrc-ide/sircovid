context("carehomes")

test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)

  ## Regenerate with: dput_named_matrix(res)
  expected <-
    rbind(icu                            = c(3, 4, 5, 10, 3),
          general                        = c(35, 34, 35, 41, 7),
          deaths_comm                    = c(23349, 23597, 23357, 23379,
                                             23245),
          deaths_hosp                    = c(283383, 283433, 282859,
                                             282832, 284127),
          admitted                       = c(131981, 132378, 132257,
                                             131758, 132269),
          new                            = c(420336, 421268, 419989,
                                             420711, 421502),
          sero_pos                       = c(3043375, 3303884, 3660340,
                                             3555769, 2384929),
          sympt_cases                    = c(12977835, 12976953, 12983125,
                                             12974885, 12981657),
          sympt_cases_over25             = c(10087409, 10086355, 10092710,
                                             10087814, 10091423),
          sympt_cases_non_variant_over25 = c(10087409, 10086355, 10092710,
                                             10087814, 10091423),
          react_pos                      = c(470, 641, 1063, 877, 151))
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

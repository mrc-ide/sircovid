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
    rbind(icu                            = c(8, 27, 7, 9, 1),
          general                        = c(30, 131, 23, 42, 10),
          deaths_comm                    = c(23673, 23552, 23356, 23393,
                                             23559),
          deaths_hosp                    = c(283863, 283386, 282999,
                                             283381, 282816),
          admitted                       = c(132186, 132127, 132194,
                                             132371, 131876),
          new                            = c(420641, 420342, 421052,
                                             420902, 420616),
          sero_pos                       = c(3479733, 4526778, 3009279,
                                             3530181, 2305953),
          sympt_cases                    = c(12981502, 12975660, 12980511,
                                             12978842, 12976269),
          sympt_cases_over25             = c(10089657, 10087311, 10087309,
                                             10089215, 10085237),
          sympt_cases_non_variant_over25 = c(10089657, 10087311, 10087309,
                                             10089215, 10085237),
          react_pos                      = c(819, 2512, 436, 836, 170))
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

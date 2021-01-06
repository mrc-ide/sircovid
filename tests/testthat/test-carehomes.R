context("carehomes")

test_that("can run the carehomes model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  end <- sircovid_date("2020-07-31") / p$dt

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  mod$set_index(carehomes_index(mod$info())$run)
  res <- mod$run(end)
saveRDS(res,"res.rds")
  expected <- rbind(
    icu = c(5, 1, 7, 4, 1),
    general = c(19, 16, 33, 32, 9),
    deaths_comm = c(23689, 23623, 23556, 23378, 23709),
    deaths_hosp = c(282876, 283382, 281952, 282694, 283949),
    admitted = c(131773, 132513, 131915, 132682, 132109),
    new = c(421294, 421505, 420611, 419953, 420487),
    sero_pos = c(3048434, 2730483, 3404946, 3360121, 2324274),
    sympt_cases = c(12977477, 12976868, 12972694, 12973542, 12972492),
    sympt_cases_over25 = c(10086654, 10087959, 10083049, 10083878, 10083429),
    react_pos = c(500, 292, 763, 685, 180))
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

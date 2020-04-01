context("parameters")

test_that("sampler runs without error", {
  time_steps_per_day <- 4
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  d <- particle_filter_data(data, "2020-02-02", time_steps_per_day)

  pars_model <- generate_parameters(
    transmission_model = "POLYMOD",
    beta = rep(0.125, 3),
    age_limits = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80) - 1,
    infection_seeding = c(0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    #severity_data_file = "extdata/severity.csv", ## the format of this file is rather different
    progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
    gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5),
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/time_steps_per_day
  )

  pars_model$beta_list <- NULL
  pars_model$beta_dates <- NULL  
  
  cmp <- readRDS("reference_pars.rds")
  expect_mapequal(pars_model, cmp)
})

context("parameters")

test_that("Parameters are generated as before", {
  time_steps_per_day <- 4

  pars_model <- generate_parameters(
    transmission_model = "POLYMOD",
    beta = rep(0.125, 3),
    progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
    gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5),
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/time_steps_per_day,
    use_polymod_pop = TRUE
  )

  # Remove newly generated parameters
  pars_model$beta_t <- NULL
  pars_model$beta_y <- NULL
  pars_model$age_bins <- NULL

  cmp <- readRDS("reference_pars.rds")

  # Remove parameters no longer generated
  pars_model$beta <- NULL
  cmp$beta <- NULL

  expect_mapequal(pars_model, cmp)
})

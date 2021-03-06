context("carehomes (support)")

test_that("carehomes progression parameters", {
  p <- carehomes_parameters_progression(0.25)
  expect_setequal(
    names(p),
    c("k_E", "k_A", "k_P", "k_C_1", "k_C_2", "k_G_D", "k_H_D", "k_H_R",
      "k_ICU_D", "k_ICU_W_R", "k_ICU_W_D", "k_ICU_pre", "k_W_D",
      "k_W_R", "k_sero_pre_1", "k_sero_pre_2", "k_sero_pos_1", "k_sero_pos_2",
      "k_PCR_pos", "k_PCR_pre", "gamma_E_step", "gamma_A_step", "gamma_P_step",
      "gamma_C_1_step", "gamma_C_2_step", "gamma_G_D_step", "gamma_H_D_step",
      "gamma_H_R_step", "gamma_ICU_D_step", "gamma_ICU_W_R_step",
      "gamma_ICU_W_D_step", "gamma_ICU_pre_step", "gamma_W_D_step",
      "gamma_W_R_step", "gamma_sero_pos_1", "gamma_sero_pos_2",
      "gamma_sero_pre_1", "gamma_sero_pre_2", "gamma_U", "gamma_PCR_pos",
      "gamma_PCR_pre", "n_gamma_E_steps", "n_gamma_A_steps", "n_gamma_P_steps",
      "n_gamma_C_1_steps", "n_gamma_C_2_steps", "n_gamma_H_D_steps",
      "n_gamma_H_R_steps", "n_gamma_ICU_D_steps", "n_gamma_ICU_W_R_steps",
      "n_gamma_ICU_W_D_steps", "n_gamma_ICU_pre_steps", "n_gamma_W_D_steps",
      "n_gamma_W_R_steps", "n_gamma_G_D_steps"))

  ## TODO: Lilith; you had said that there were some constraints
  ## evident in the fractional representation of these values - can
  ## you write tests here that reflect that?
})


test_that("carehomes vaccination parameters", {
  n_groups <- carehomes_n_groups()
  # test default values
  ntot <- rep(1000, n_groups)
  p <- carehomes_parameters_vaccination(ntot)
  expect_setequal(
    names(p),
    c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_death",
      "rel_infectivity",
      "n_vacc_classes",
      "vaccine_progression_rate_base", "vaccine_dose_step",
      "vaccine_catchup_fraction",
      "index_dose", "index_dose_inverse", "n_doses"))
  expect_equal(nrow(p$rel_susceptibility), n_groups)
  expect_equal(ncol(p$rel_susceptibility), 1)
  expect_equal(nlayer(p$rel_susceptibility), 1)
  expect_equal(nrow(p$vaccine_progression_rate_base), n_groups)
  expect_equal(ncol(p$vaccine_progression_rate_base), 1)

  # test when more vaccinated categories than default
  rel_susceptibility <- c(1, 0.2, 0.1, 0.4)
  rel_p_sympt <- c(1, 0.75, 0.5, 0.75)
  rel_p_hosp_if_sympt <- c(1, 0.8, 0.6, 0.9)
  vaccine_progression_rate <- c(0, 1, 0, 1)

  region <- "london"
  vaccine_daily_doses <- c(5000, 10000)
  vaccine_daily_doses_date <- sircovid_date(c("2020-12-01", "2021-02-01"))
  daily_doses <- c(rep(vaccine_daily_doses[1], diff(vaccine_daily_doses_date)),
                   rep(vaccine_daily_doses[2], 200))
  uptake <- test_example_uptake()
  n <- vaccine_priority_population(region, uptake)

  vaccine_schedule <- vaccine_schedule_future(
    vaccine_daily_doses_date[1], daily_doses, 1e6, n)

  p <- carehomes_parameters_vaccination(ntot,
                                        dt = 0.25,
                                        rel_susceptibility = rel_susceptibility,
                                        rel_p_sympt = rel_p_sympt,
                                        rel_p_hosp_if_sympt =
                                          rel_p_hosp_if_sympt,
                                        vaccine_progression_rate =
                                          vaccine_progression_rate,
                                        vaccine_schedule = vaccine_schedule,
                                        vaccine_index_dose2 = 3L)
  expect_setequal(
    names(p),
    c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_death",
      "rel_infectivity",
      "n_vacc_classes",
      "vaccine_progression_rate_base", "vaccine_dose_step",
      "vaccine_catchup_fraction",
      "index_dose", "index_dose_inverse", "n_doses"))

  expect_equal(nrow(p$rel_susceptibility), n_groups)
  expect_equal(ncol(p$rel_susceptibility), 1)
  expect_equal(nlayer(p$rel_susceptibility), length(rel_susceptibility))

  expect_equal(nrow(p$rel_p_sympt), n_groups)
  expect_equal(ncol(p$rel_p_sympt), 1)
  expect_equal(dim(p$rel_p_sympt)[3], length(rel_p_sympt))

  expect_equal(nrow(p$rel_p_hosp_if_sympt), n_groups)
  expect_equal(ncol(p$rel_p_hosp_if_sympt), 1)
  expect_equal(dim(p$rel_p_hosp_if_sympt)[3], length(rel_p_hosp_if_sympt))

  expect_equal(nrow(p$rel_infectivity), n_groups)
  expect_equal(ncol(p$rel_infectivity), 1)
  expect_equal(dim(p$rel_infectivity)[3], length(rel_susceptibility))

  expect_equal(nrow(p$vaccine_progression_rate_base), n_groups)
  expect_equal(ncol(p$vaccine_progression_rate_base),
               length(vaccine_progression_rate))

  expect_equal(dim(p$vaccine_dose_step),
               c(19, 2,
                 (length(daily_doses) + vaccine_daily_doses_date[1]) * 4))
  ## daily doses are as expected
  expect_vector_equal(colSums(p$vaccine_dose_step[, 1, ]),
                      c(rep(0, vaccine_daily_doses_date[1] * 4),
                        rep(daily_doses / 4, each = 4)),
                      digits = 0, tol = 2)
  msg1 <- "rel_susceptibility, rel_p_sympt, rel_p_hosp_if_sympt, rel_p_death,"
  msg2 <- "rel_infectivity should have the same dimension"
  expect_error(
    carehomes_parameters_vaccination(ntot,
                                     rel_susceptibility = 1,
                                     rel_p_sympt = c(1, 0.5, 0.25),
                                     rel_p_hosp_if_sympt = c(1, 0.1),
                                     rel_p_death = c(1, 0.2),
                                     rel_infectivity = 1),
    paste(msg1, msg2))
  expect_error(carehomes_parameters_vaccination(ntot,
                                                rel_susceptibility = c(1, 1),
                                                rel_p_sympt = c(1, 0.5, 0.25),
                                                rel_p_hosp_if_sympt = 1,
                                                rel_p_death = 1,
                                                rel_infectivity = 1),
               paste(msg1, msg2))
  expect_error(carehomes_parameters_vaccination(ntot,
                                                rel_susceptibility = c(1, 1),
                                                rel_p_sympt = c(1, 0.5),
                                                rel_p_hosp_if_sympt = c(1, 1),
                                                rel_p_death = c(1, 0.2),
                                                rel_infectivity = c(1, 1, 0.5)),
               paste(msg1, msg2))
  expect_error(carehomes_parameters_vaccination(ntot,
                                                vaccine_catchup_fraction = -1),
               "'vaccine_catchup_fraction' must lie in [0, 1]",
               fixed = TRUE)
  expect_error(carehomes_parameters_vaccination(ntot,
                                                vaccine_catchup_fraction = 1.5),
               "'vaccine_catchup_fraction' must lie in [0, 1]",
               fixed = TRUE)
  expect_error(
    carehomes_parameters_vaccination(ntot,
                                     vaccine_catchup_fraction = c(1, 1)),
    "'vaccine_catchup_fraction' must be a scalar")
})

test_that("carehomes_parameters returns a list of parameters", {
  date <- sircovid_date("2020-02-01")
  p <- carehomes_parameters(date, "uk")

  expect_type(p, "list")
  expect_equal(p$beta_step, 0.08)

  ## Transmission matrix is more complicated - see below
  expect_identical(p$m[1:17, 1:17], sircovid_transmission_matrix("uk"))
  expect_identical(
    p$m,
    carehomes_transmission_matrix(0.1, 4e-6, 5e-5, "uk"))

  progression <- carehomes_parameters_progression(0.25)
  expect_identical(p[names(progression)], progression)

  vaccination <- carehomes_parameters_vaccination(
    p$N_tot, p$dt, p$rel_susceptibility, p$rel_p_sympt, p$rel_p_hosp_if_sympt,
    p$rel_p_death, p$rel_infectivity,
    p$vaccine_progression_rate_base)
  expect_identical(p[names(vaccination)], vaccination)

  strain <- carehomes_parameters_strain(p$strain_transmission,
                                        strain_seed_date = NULL,
                                        strain_seed_rate = NULL,
                                        dt = 1 / 4)

  waning <- carehomes_parameters_waning(0)
  expect_identical(p[names(waning)], waning)

  shared <- sircovid_parameters_shared(date, "uk", NULL, NULL)
  expect_identical(p[names(shared)], shared)

  severity <- carehomes_parameters_severity(0.25, NULL)
  expect_identical(p[names(severity)], severity)

  observation <- carehomes_parameters_observation(1e6)
  expect_identical(p[names(observation)], observation)

  sens_and_spec <- carehomes_parameters_sens_and_spec()
  expect_identical(p[names(sens_and_spec)], sens_and_spec)

  expect_equal(p$N_tot_15_64, sum(p$N_tot[4:13]))

  extra <- setdiff(names(p),
                   c("m", names(observation),
                     names(shared), names(progression), names(severity),
                     names(strain), names(vaccination), names(waning),
                     names(sens_and_spec)))
  expect_setequal(
    extra,
    c("N_tot", "carehome_beds", "carehome_residents", "carehome_workers",
      "rel_p_ICU", "rel_p_ICU_D", "rel_p_H_D", "rel_p_W_D", "rel_p_G_D",
      "rel_p_R", "rel_gamma_E", "rel_gamma_A", "rel_gamma_P", "rel_gamma_C_1",
      "rel_gamma_C_2", "rel_gamma_H_D", "rel_gamma_H_R", "rel_gamma_ICU_pre",
      "rel_gamma_ICU_D", "rel_gamma_ICU_W_D", "rel_gamma_ICU_W_R",
      "rel_gamma_W_D", "rel_gamma_W_R", "rel_gamma_G_D",
      "N_tot_15_64", "N_tot_all", "N_tot_over25", "N_tot_react",
      "p_NC", "I_A_transmission", "I_P_transmission",
      "I_C_1_transmission", "I_C_2_transmission",
      "n_groups", "initial_I", "cross_immunity"))

  expect_equal(p$carehome_beds, sircovid_carehome_beds("uk"))
  expect_equal(p$carehome_residents, round(p$carehome_beds * 0.742))
  expect_equal(p$carehome_workers, p$carehome_residents)

  ## Can be slightly off due to rounding error
  expect_true(abs(sum(p$N_tot) - sum(p$population)) < 3)
  expect_length(p$N_tot, 19)
})


test_that("can compute severity for carehomes model", {
  severity <- carehomes_parameters_severity(0.25, NULL)

  expect_equal(
    severity$p_G_D_step, array(0.05, c(1, 19)))
  expect_equal(
    severity$p_star_step, array(0.2, c(1, 19)))
})


test_that("can input severity data for carehomes model", {
  data <- sircovid_parameters_severity(NULL)

  data$p_G_D[] <- 0
  data$p_star[] <- 0.5

  severity <- carehomes_parameters_severity(0.25, data)

  expect_equal(
    severity$p_G_D_step, array(0, c(1, 19)))
  expect_equal(
    severity$p_star_step, array(0.5, c(1, 19)))
})


test_that("can compute time-varying severity parameters for carehomes model", {
  dt <- 0.25

  p_G_D_date <- sircovid_date(c("2020-02-01", "2020-05-01"))
  p_G_D_value <- c(0.05, 0.1)
  p_G_D_CHR_value <- 0.4

  p_H_value <- 0.6
  p_H_CHR_date <- sircovid_date(c("2020-03-01", "2020-04-01"))
  p_H_CHR_value <- c(0.7, 0.6)

  severity <-
    carehomes_parameters_severity(dt, NULL,
                                  p_H = list(value = p_H_value),
                                  p_H_CHR = list(date = p_H_CHR_date,
                                                 value = p_H_CHR_value),
                                  p_G_D = list(date = p_G_D_date,
                                               value = p_G_D_value),
                                  p_G_D_CHR = list(value = p_G_D_CHR_value))

  p_G_D_step <-
    sircovid_parameters_piecewise_linear(p_G_D_date,
                                         p_G_D_value, dt)
  expect_equal(severity$p_G_D_step[, 19],
               rep(p_G_D_CHR_value, length(p_G_D_step)))
  expect_equal(severity$p_G_D_step[, 17], p_G_D_step)


  p_H_CHR_step <-
    sircovid_parameters_piecewise_linear(p_H_CHR_date,
                                         p_H_CHR_value, dt)
  expect_equal(severity$p_H_step[, 17],
               rep(p_H_value, length(p_H_CHR_step)))
  expect_equal(severity$p_H_step[, 19], p_H_CHR_step)

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_C = list(date = 1,
                                             value = 0.3)),
    "As 'p_C' has a single 'value', expected NULL or missing 'date'")

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_ICU = list(date = c(1, 4, 5),
                                               value = c(0.2, 0.3))),
    "'date' and 'value' for 'p_ICU' must have the same length")

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_ICU_D = list(date = c(1, 4),
                                                 value = c(-1, 0.3))),
    "'p_ICU_D' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_W_D = list(date = c(1, 4),
                                               value = c(0.2, 3))),
    "'p_W_D' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_H_CHR = list(date = 1,
                                                 value = 0.3)),
    "As 'p_H_CHR' has a single 'value', expected NULL or missing 'date'")

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_G_D_CHR = list(date = c(1, 4, 5),
                                                   value = c(0.2, 0.3))),
    "'date' and 'value' for 'p_G_D_CHR' must have the same length")

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_H_CHR = list(date = c(1, 4),
                                                 value = c(-1, 0.3))),
    "'p_H_CHR' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    carehomes_parameters_severity(dt,
                                  p_G_D_CHR = list(date = c(1, 4),
                                                   value = c(0.2, 3))),
    "'p_G_D_CHR' must lie in [0, 1]", fixed = TRUE)

})


test_that("can compute time-varying progression parameters for carehomes
          model", {
  dt <- 0.25

  gamma_H_R_value <- 0.3
  gamma_H_D_date <- sircovid_date(c("2020-02-01", "2020-05-01"))
  gamma_H_D_value <- c(0.2, 0.5)

  progression <-
    carehomes_parameters_progression(dt,
                                     gamma_H_D = list(date = gamma_H_D_date,
                                                      value = gamma_H_D_value),
                                     gamma_H_R = list(value = gamma_H_R_value)
    )

  gamma_H_D_step <-
    sircovid_parameters_piecewise_linear(gamma_H_D_date,
                                         gamma_H_D_value, dt)
  expect_equal(progression$gamma_H_D_step, gamma_H_D_step)
  expect_equal(progression$gamma_H_R_step, gamma_H_R_value)

  expect_error(
    carehomes_parameters_progression(dt,
                                     gamma_E = list(date = 1,
                                                    value = 3)),
    "As 'gamma_E' has a single 'value', expected NULL or missing 'date'")

  expect_error(
    carehomes_parameters_progression(dt,
                                     gamma_ICU_pre = list(date = c(1, 4, 5),
                                                          value = c(2, 3))),
    "'date' and 'value' for 'gamma_ICU_pre' must have the same length")

  expect_error(
    carehomes_parameters_progression(dt,
                                     gamma_H_D = list(date = c(1, 4),
                                                      value = c(-2, 3))),
    "'gamma_H_D' must have only non-negative values")


})


test_that("Can compute transmission matrix for carehomes model", {
  m <- carehomes_transmission_matrix(0.1, 4e-5, 5e-4, "uk")
  expect_equal(rownames(m)[18:19], c("CHW", "CHR"))
  expect_equal(colnames(m)[18:19], c("CHW", "CHR"))
  expect_equal(dim(m), c(19, 19))
  expect_equal(m[1:17, 1:17], sircovid_transmission_matrix("uk"))
  expect_equal(unname(m[18:19, 18:19]),
               matrix(rep(c(4e-5, 5e-4), c(3, 1)), 2))
  expect_equal(m, t(m))

  ## TODO: we should get a check of the weighted mean and resident
  ## scaling here
})


test_that("can tune the noise parameter", {
  p1 <- carehomes_parameters_observation(1e6)
  p2 <- carehomes_parameters_observation(1e4)
  expect_setequal(names(p1), names(p2))
  v <- setdiff(names(p1), "exp_noise")
  expect_mapequal(p1[v], p2[v])
  expect_equal(p1$exp_noise, 1e6)
  expect_equal(p2$exp_noise, 1e4)
})


test_that("carehomes_index identifies ICU and D_tot in real model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 10)
  info <- mod$info()
  index <- carehomes_index(info)

  expect_equal(
    names(index$run),
    c("icu", "general", "deaths_carehomes_inc", "deaths_comm_inc",
      "deaths_hosp_inc", "admitted_inc", "diagnoses_inc",
      "sero_pos_1", "sero_pos_2", "sympt_cases_inc",
      "sympt_cases_non_variant_inc", "sympt_cases_over25_inc",
      "sympt_cases_non_variant_over25_inc", "react_pos"))

  expect_equal(index$run[["icu"]],
               which(names(info$index) == "ICU_tot"))
  expect_equal(index$run[["general"]],
               which(names(info$index) == "general_tot"))
  expect_equal(index$run[["deaths_carehomes_inc"]],
               which(names(info$index) == "D_carehomes_inc"))
  expect_equal(index$run[["deaths_comm_inc"]],
               which(names(info$index) == "D_comm_inc"))
  expect_equal(index$run[["deaths_hosp_inc"]],
               which(names(info$index) == "D_hosp_inc"))
  expect_equal(index$run[["admitted_inc"]],
               which(names(info$index) == "admit_conf_inc"))
  expect_equal(index$run[["diagnoses_inc"]],
               which(names(info$index) == "new_conf_inc"))
  expect_equal(index$run[["sero_pos_1"]],
               which(names(info$index) == "sero_pos_1"))
  expect_equal(index$run[["sero_pos_2"]],
               which(names(info$index) == "sero_pos_2"))
  expect_equal(index$run[["sympt_cases_inc"]],
               which(names(info$index) == "sympt_cases_inc"))
  expect_equal(index$run[["sympt_cases_over25_inc"]],
               which(names(info$index) == "sympt_cases_over25_inc"))
  expect_equal(index$run[["react_pos"]],
               which(names(info$index) == "react_pos"))
})


test_that("carehome worker index is sensible", {
  expect_equal(carehomes_index_workers(), 6:13) })


test_that("carehomes_index is properly named", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 10)
  info <- mod$info()
  index <- carehomes_index(info)
  expect_false(any(is.na(names(index))))
})


test_that("Can compute initial conditions", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 10)
  info <- mod$info()

  initial <- carehomes_initial(info, 10, p)
  expect_setequal(names(initial), c("state", "step"))
  expect_equal(initial$step, p$initial_step)

  initial_y <- mod$transform_variables(initial$state)

  expect_equal(initial_y$N_tot_sero_1, sum(p$N_tot))
  expect_equal(initial_y$N_tot_sero_2, sum(p$N_tot))
  expect_equal(initial_y$N_tot_PCR, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S) + drop(initial_y$I_A),
               p$N_tot)
  expect_equal(drop(initial_y$I_A),
               append(rep(0, 18), 10, after = 3))
  expect_equal(drop(initial_y$I_weighted),
               append(rep(0, 18), p$I_A_transmission * 10, after = 3))
  expect_equal(initial_y$T_sero_pre_1[, 1, 1, ],
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$T_sero_pre_2[, 1, 1, ],
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$T_PCR_pos[, 1, 1, ],
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$react_pos, 10)

  ## 48 here, derived from;
  ## * 38 (S + N_tot)
  ## * 1 (prob_strain)
  ## * 1 (react_pos)
  ## * 3 (N_tot_sero_1 + N_tot_sero_2 + N_tot_PCR)
  ## * 2 (I_A[4] + I_weighted[4])
  ## * 3 (T_sero_pre_1[4] + T_sero_pre_2[4] + T_PCR_pos[4])
  expect_equal(sum(initial$state != 0), 48)
})


test_that("Can control the seeding", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "london",
                            initial_I = 50)
  expect_equal(p$initial_I, 50)

  mod <- carehomes$new(p, 0, 10)
  info <- mod$info()

  initial <- carehomes_initial(info, 10, p)
  expect_setequal(names(initial), c("state", "step"))
  expect_equal(initial$step, p$initial_step)

  initial_y <- mod$transform_variables(initial$state)

  expect_equal(initial_y$N_tot_sero_1, sum(p$N_tot))
  expect_equal(initial_y$N_tot_sero_2, sum(p$N_tot))
  expect_equal(initial_y$N_tot_PCR, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S) + drop(initial_y$I_A),
               p$N_tot)
  expect_equal(drop(initial_y$I_A),
               append(rep(0, 18), 50, after = 3))
  expect_equal(initial_y$T_sero_pre_1[, 1, 1, ],
               append(rep(0, 18), 50, after = 3))
  expect_equal(initial_y$T_sero_pre_2[, 1, 1, ],
               append(rep(0, 18), 50, after = 3))
  expect_equal(initial_y$T_PCR_pos[, 1, 1, ],
               append(rep(0, 18), 50, after = 3))
  expect_equal(initial_y$react_pos, 50)

  ## 48 here, derived from;
  ## * 38 (S + N_tot)
  ## * 1 (prob_strain)
  ## * 1 (react_pos)
  ## * 3 (N_tot_sero_1 + N_tot_sero_2 + N_tot_PCR)
  ## * 2 (I_A[4] + I_weighted[4])
  ## * 3 (T_sero_pre_1[4] + T_sero_pre_2[4] + T_PCR_pos[4])
  expect_equal(sum(initial$state != 0), 48)
})


test_that("get carehomes population can get the population", {
  expect_equal(sircovid_carehome_beds("uk"), 537521)
  expect_equal(sircovid_carehome_beds("UK"), 537521)
})


test_that("sircovid_carehome_beds throws sensible error on invalid input", {
  expect_error(sircovid_carehome_beds(NULL), "'region' must not be NULL")
  expect_error(
    sircovid_carehome_beds("oxfordshire"),
    "Carehome beds not found for 'oxfordshire': must be one of 'east_of_")
})


test_that("sircovid_carehome_beds caches data", {
  clear_cache()
  n <- sircovid_carehome_beds("uk")
  expect_s3_class(cache$carehomes, "data.frame")
  expect_equal(cache$carehomes$carehome_beds[cache$carehomes$region == "uk"],
               n)
})


## TODO: Ed - you had said that you had ideas for some more systematic
## testing here.  This will also get easier to do if we move to having
## a function generator given a data set.
test_that("carehomes_compare combines likelihood correctly", {
  state <- rbind(
    icu = 10:15,
    general = 20:25,
    deaths_carehomes_inc = 2:7,
    deaths_comm_inc = 1:6,
    deaths_hosp_inc = 3:8,
    admitted_inc = 50:55,
    diagnoses_inc = 60:65,
    sero_pos_1 = 4:9,
    sero_pos_2 = 14:19,
    sympt_cases_inc = 100:105,
    sympt_cases_non_variant_inc = 70:75,
    sympt_cases_over25_inc = 80:85,
    sympt_cases_non_variant_over25_inc = 60:65,
    react_pos = 2:7)
  observed <- list(
    icu = 13,
    general = 23,
    hosp = 36,
    deaths_carehomes = 4,
    deaths_hosp = 5,
    deaths_comm = 3,
    deaths = 8,
    deaths_non_hosp = 6,
    admitted = 53,
    diagnoses = 63,
    all_admission = 116,
    sero_pos_15_64_1 = 43,
    sero_tot_15_64_1 = 83,
    sero_pos_15_64_2 = 58,
    sero_tot_15_64_2 = 98,
    pillar2_pos = 35,
    pillar2_tot = 600,
    pillar2_cases = 35,
    pillar2_over25_pos = 25,
    pillar2_over25_tot = 500,
    pillar2_over25_cases = 25,
    react_pos = 3,
    react_tot = 500,
    strain_non_variant = 40,
    strain_tot = 50,
    strain_over25_non_variant = 20,
    strain_over25_tot = 25)
  date <- sircovid_date("2020-01-01")
  pars <- carehomes_parameters(date, "uk", exp_noise = Inf)

  observed_keep <- function(nms) {
    observed[setdiff(names(observed), nms)] <- NA_real_
    observed
  }
  observed_drop <- function(nms) {
    observed[nms] <- NA_real_
    observed
  }

  ## This function is more complicated to test than the basic model
  ## because it's not a simple sum
  nms_sero_1 <- c("sero_pos_15_64_1", "sero_tot_15_64_1")
  nms_sero_2 <- c("sero_pos_15_64_2", "sero_tot_15_64_2")
  nms_pillar2 <- c("pillar2_pos", "pillar2_tot")
  nms_pillar2_over25 <- c("pillar2_over25_pos", "pillar2_over25_tot")
  nms_react <- c("react_pos", "react_tot")
  nms_strain <- c("strain_non_variant", "strain_tot")
  nms_strain_over25 <- c("strain_over25_non_variant", "strain_over25_tot")
  parts <- c(as.list(setdiff(names(observed),
                             c(nms_sero_1, nms_sero_2, nms_pillar2,
                               nms_pillar2_over25, nms_react, nms_strain,
                               nms_strain_over25))),
             list(nms_sero_1), list(nms_sero_2), list(nms_pillar2),
             list(nms_pillar2_over25), list(nms_react), list(nms_strain),
             list(nms_strain_over25))

  ll_parts <- lapply(parts, function(x)
    carehomes_compare(state, observed_keep(x), pars))

  ## Extremely light testing, though this has already flushed out some
  ## issues
  expect_vector_equal(lengths(ll_parts), 6)
  expect_equal(
    carehomes_compare(state, observed, pars),
    rowSums(do.call(cbind, ll_parts)))

  ## Test that there are non-zero values for each log-likelihood part.
  ## This helps make sure all parts contribute to the log-likelihood.
  expect_true(all(sapply(ll_parts, function(x) any(x != 0))))
})


test_that("carehomes_population prevents negative populations", {
  population <- rep(100, 17)
  expect_error(
    carehomes_population(population, 1000, 0),
    "Not enough population to be care workers")
  expect_error(
    carehomes_population(population, 0, 1000),
    "Not enough population to meet care home occupancy")
})


test_that("carehomes_population preserves population", {
  population <- c(949, 989, 1030, 995, 993, 985, 965, 1082, 1042, 960,
                  980, 1004, 934, 1049, 971, 1020, 937)
  res <- carehomes_population(population, 120, 200)
  expect_equal(sum(res), sum(population))
  expect_equal(res[18:19], c(120, 200))
  expect_equal(res[1:5], population[1:5])
  expect_vector_lt(res[6:17], population[6:17])
})


test_that("carehomes_index returns S compartments", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  index <- carehomes_index(mod$info())
  info <- mod$info()
  i <- seq_along(index$run)
  expect_equal(
    unname(index$state[!names(index$state) %in% names(index$run)]),
    c(info$index$D_comm_tot,
      info$index$D_carehomes_tot,
      info$index$D_hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$cum_sympt_cases,
      info$index$cum_sympt_cases_non_variant,
      info$index$cum_sympt_cases_over25,
      info$index$cum_sympt_cases_non_variant_over25,
      info$index$hosp_tot,
      info$index$D_tot,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$S,
      info$index$R,
      info$index$cum_admit_by_age,
      info$index$D_hosp,
      info$index$D,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$I_weighted,
      info$index$prob_strain,
      info$index$cum_n_vaccinated))
})


test_that("carehomes_particle_filter_data requires consistent deaths", {
  data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                        1, 0.25)
  ## Add additional columns
  data$deaths_hosp <- data$deaths
  data$deaths_carehomes <- NA
  data$deaths_comm <- NA
  data$deaths_non_hosp <- NA
  data$general <- NA
  data$hosp <- NA
  data$admitted <- NA
  data$diagnoses <- NA
  data$all_admission <- NA
  data$sero_pos_15_64_1 <- NA
  data$sero_tot_15_64_1 <- NA
  data$sero_pos_15_64_2 <- NA
  data$sero_tot_15_64_2 <- NA
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
  data$strain_over25_non_variant <- NA
  data$strain_over25_tot <- NA
  expect_error(
    carehomes_particle_filter(data),
    "Deaths are not consistently split into total vs hospital/non-hospital
          or hospital/care homes/community")

  data$deaths_non_hosp <- data$deaths
  data$deaths_carehomes <- data$deaths
  data$deaths <- NA
  expect_error(
    carehomes_particle_filter(data),
    "Non-hospital deaths are not consistently split into total vs care
         homes/community")
})


test_that("carehomes_particle_filter_data does not allow more than one pillar 2
          or strain data stream", {
            data <- sircovid_data(
              read_csv(sircovid_file("extdata/example.csv")), 1, 0.25)
            ## Add additional columns
            data$deaths_hosp <- NA
            data$deaths_comm <- NA
            data$deaths_carehomes <- NA
            data$deaths_non_hosp <- NA
            data$general <- NA
            data$hosp <- NA
            data$admitted <- NA
            data$diagnoses <- NA
            data$all_admission <- NA
            data$sero_pos_15_64_1 <- NA
            data$sero_tot_15_64_1 <- NA
            data$sero_pos_15_64_2 <- NA
            data$sero_tot_15_64_2 <- NA
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
            data$strain_over25_non_variant <- NA
            data$strain_over25_tot <- NA

            ## Example that should not cause error
            data$pillar2_pos <- 3
            data$pillar2_tot <- 6
            data$strain_non_variant <- 8
            data$strain_tot <- 10
            pf <- carehomes_particle_filter(data, 10)
            expect_s3_class(pf, "particle_filter")

            data$pillar2_pos <- 3
            data$pillar2_tot <- 6
            data$pillar2_cases <- 5
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one pillar 2 data stream")

            data$pillar2_pos <- 3
            data$pillar2_tot <- 6
            data$pillar2_cases <- NA
            data$pillar2_over25_cases <- 5
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one pillar 2 data stream")

            data$pillar2_pos <- 3
            data$pillar2_tot <- 6
            data$pillar2_over25_cases <- 5
            data$pillar2_over25_pos <- 3
            data$pillar2_over25_tot <- 6
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one pillar 2 data stream")

            data$pillar2_pos <- NA
            data$pillar2_tot <- NA
            data$pillar2_over25_cases <- 5
            data$pillar2_over25_pos <- 3
            data$pillar2_over25_tot <- 6
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one pillar 2 data stream")

            data$pillar2_over25_cases <- NA
            data$strain_non_variant <- 8
            data$strain_tot <- 10
            data$strain_over25_non_variant <- 6
            data$strain_over25_tot <- 9
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one strain data stream")

          })


test_that("the carehomes sircovid model has 19 groups", {
  expect_equal(carehomes_n_groups(), 19)
})


test_that("carehomes_particle_filter_data does not allow more than one pillar 2
          or strain data stream (II)", {
            data <- sircovid_data(
              read_csv(sircovid_file("extdata/example.csv")), 1, 0.25)
            class(data)[1] <- "particle_filter_data_nested"
            ## Add additional columns
            data$deaths_hosp <- NA
            data$deaths_comm <- NA
            data$deaths_carehomes <- NA
            data$deaths_non_hosp <- NA
            data$general <- NA
            data$hosp <- NA
            data$admitted <- NA
            data$diagnoses <- NA
            data$all_admission <- NA
            data$sero_pos_15_64_1 <- NA
            data$sero_tot_15_64_1 <- NA
            data$sero_pos_15_64_2 <- NA
            data$sero_tot_15_64_2 <- NA
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
            data$strain_over25_non_variant <- NA
            data$strain_over25_tot <- NA

            ## Example that should not cause error
            data$pillar2_pos <- 3
            data$pillar2_tot <- 6

            ## Add populations
            data <- rbind(data, data)
            data$population <- rep(factor(letters[1:2]), each = 31)
            pf <- carehomes_particle_filter(data, 10)
            expect_s3_class(pf, "particle_filter")

            ## Error if one region has multiple streams
            data$pillar2_pos <- 3
            data$pillar2_tot <- 6
            data$pillar2_cases <- 5
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one pillar 2 data stream")

            data$pillar2_cases <- NA
            data$strain_non_variant <- 8
            data$strain_tot <- 10
            data$strain_over25_non_variant <- 6
            data$strain_over25_tot <- 9
            expect_error(
              carehomes_particle_filter(data),
              "Cannot fit to more than one strain data stream")

            ## Don't error if two regions have different streams
            data[1:31, "pillar2_cases"] <- NA
            data[32:62, "pillar2_pos"] <- NA
            data[32:62, "pillar2_tot"] <- NA
            data[1:31, "strain_non_variant"] <- NA
            data[1:31, "strain_tot"] <- NA
            data[32:62, "strain_over25_non_variant"] <- NA
            data[32:62, "strain_over25_tot"] <- NA
            pf <- carehomes_particle_filter(data, 10)
            expect_s3_class(pf, "particle_filter")
          })


test_that("carehomes check severity works as expected", {
  ## errors if params missing
  expect_error(
    carehomes_check_severity(list(n_groups = 19)),
    "Parameter 'rel_p_sympt' is missing"
  )
  expect_error(
    carehomes_check_severity(list(n_groups = 19, rel_p_sympt = 1)),
    "Parameter 'p_C_step' is missing"
  )

  ## errors if rel_ is not 1 or 4 cols
  expect_error(carehomes_check_severity(list(
    n_groups = 19,
    rel_p_sympt = array(1, c(19, 3, 3)),
    p_C_step = matrix(1, 1, 19)
  )), "1 or 4 columns")

  ## no error if 1 col
  expect_error(carehomes_check_severity(list(
    n_groups = 19,
    rel_p_sympt = array(runif(19 * 3), c(19, 1, 3)),
    p_C_step = matrix(1, 1, 19)
  )), "Parameter 'rel_p_hosp_if_sympt' is missing", fixed = TRUE)

  ## check on required parameters
  steps <- c(
    "p_C_step", "p_H_step", "p_ICU_step", "p_ICU_D_step",
    "p_H_D_step", "p_W_D_step", "p_G_D_step", "p_R_step"
  )
  rels <- c(
    "rel_p_sympt",
    "rel_p_hosp_if_sympt", "rel_p_ICU", "rel_p_ICU_D",
    "rel_p_H_D", "rel_p_W_D", "rel_p_G_D", "rel_p_R"
  )

  p <- vector("list", 16)
  names(p) <- c(steps, rels)
  ## 19 groups, 4 strains, 3 vacc classes
  p[1:8] <- list(matrix(0.5, ncol = 19))
  p[9:16] <- list(array(1, c(19, 4, 3)))

  p$n_groups <- 19

  ## first check no errors
  expect_equal(carehomes_check_severity(p), p)

  ## check errors when expected - we only check upper bound as it's safe to
  ##  assume negative probs won't be given in practice (but uses the same )

  ## rel_p_sympt Inf
  for (i in seq_along(steps)) {
    p[[rels[[i]]]][] <- Inf
    expect_error(
      carehomes_check_severity(p),
      sprintf("%s * %s is not in [0, 1] for group 1", rels[[i]], steps[[i]]),
      fixed = TRUE
    )
    p[[rels[[i]]]][] <- 1
  }
})


test_that("Can input population data", {

  pop <- sircovid_population("uk")
  p <- carehomes_parameters(1, "Unknown", population = pop, carehome_beds = 0)
  expect_equal(p$population, pop)
  expect_equal(p$N_tot, c(pop, 0, 0))

  expect_error(
    carehomes_parameters(1, "Unknown", population = pop[-1L],
                         carehome_beds = 0),
    "If population is specified it must be a vector of length 17")

  pop[1L] <- 0.5
  expect_error(
    carehomes_parameters(1, "Unknown", population = pop, carehome_beds = 0),
    "'population' must be an integer")

  pop[1L] <- -1
  expect_error(
    carehomes_parameters(1, "Unknown", population = pop, carehome_beds = 0),
    "'population' must have only non-negative values")
})

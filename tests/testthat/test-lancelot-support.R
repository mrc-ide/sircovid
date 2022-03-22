context("lancelot (support)")

test_that("lancelot progression parameters", {
  p <- lancelot_parameters_progression(0.25)
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
      "gamma_sero_pre_1", "gamma_sero_pre_2", "gamma_U_step", "gamma_PCR_pos",
      "gamma_PCR_pre", "n_gamma_E_steps", "n_gamma_A_steps", "n_gamma_P_steps",
      "n_gamma_C_1_steps", "n_gamma_C_2_steps", "n_gamma_H_D_steps",
      "n_gamma_H_R_steps", "n_gamma_ICU_D_steps", "n_gamma_ICU_W_R_steps",
      "n_gamma_ICU_W_D_steps", "n_gamma_ICU_pre_steps", "n_gamma_W_D_steps",
      "n_gamma_W_R_steps", "n_gamma_G_D_steps", "n_gamma_U_steps"))

  ## TODO: Lilith; you had said that there were some constraints
  ## evident in the fractional representation of these values - can
  ## you write tests here that reflect that?
})


test_that("lancelot vaccination parameters", {
  n_groups <- lancelot_n_groups()
  # test default values
  ntot <- rep(1000, n_groups)
  p <- lancelot_parameters_vaccination(ntot)
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

  p <- lancelot_parameters_vaccination(ntot,
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
    lancelot_parameters_vaccination(ntot,
                                    rel_susceptibility = 1,
                                    rel_p_sympt = c(1, 0.5, 0.25),
                                    rel_p_hosp_if_sympt = c(1, 0.1),
                                    rel_p_death = c(1, 0.2),
                                    rel_infectivity = 1),
    paste(msg1, msg2))
  expect_error(lancelot_parameters_vaccination(ntot,
                                               rel_susceptibility = c(1, 1),
                                               rel_p_sympt = c(1, 0.5, 0.25),
                                               rel_p_hosp_if_sympt = 1,
                                               rel_p_death = 1,
                                               rel_infectivity = 1),
               paste(msg1, msg2))
  expect_error(lancelot_parameters_vaccination(ntot,
                                               rel_susceptibility = c(1, 1),
                                               rel_p_sympt = c(1, 0.5),
                                               rel_p_hosp_if_sympt = c(1, 1),
                                               rel_p_death = c(1, 0.2),
                                               rel_infectivity = c(1, 1, 0.5)),
               paste(msg1, msg2))
  expect_error(lancelot_parameters_vaccination(ntot,
                                               vaccine_catchup_fraction = -1),
               "'vaccine_catchup_fraction' must lie in [0, 1]",
               fixed = TRUE)
  expect_error(lancelot_parameters_vaccination(ntot,
                                               vaccine_catchup_fraction = 1.5),
               "'vaccine_catchup_fraction' must lie in [0, 1]",
               fixed = TRUE)
  expect_error(
    lancelot_parameters_vaccination(ntot,
                                    vaccine_catchup_fraction = c(1, 1)),
    "'vaccine_catchup_fraction' must be a scalar")
})

test_that("lancelot_parameters returns a list of parameters", {
  date <- sircovid_date("2020-02-01")
  p <- lancelot_parameters(date, "uk")

  expect_type(p, "list")
  expect_equal(p$beta_step, 0.1)

  ## Transmission matrix is more complicated - see below
  expect_identical(p$m[1:17, 1:17], sircovid_transmission_matrix("uk"))
  expect_identical(
    p$m,
    lancelot_transmission_matrix(0.1, 4e-6, 5e-5, "uk"))

  progression <- lancelot_parameters_progression(0.25)
  expect_identical(p[names(progression)], progression)

  vaccination <- lancelot_parameters_vaccination(
    p$N_tot, p$dt, p$rel_susceptibility, p$rel_p_sympt, p$rel_p_hosp_if_sympt,
    p$rel_p_death, p$rel_infectivity,
    p$vaccine_progression_rate_base)
  expect_identical(p[names(vaccination)], vaccination)

  strain <- lancelot_parameters_strain(p$strain_transmission,
                                       strain_seed_date = NULL,
                                       strain_seed_size = NULL,
                                       strain_seed_pattern = NULL,
                                       dt = 1 / 4)

  waning <- lancelot_parameters_waning(0)
  expect_identical(p[names(waning)], waning)

  shared <- sircovid_parameters_shared(date, "uk", NULL, NULL,
                                       "piecewise-linear", NULL, 1, 30)
  expect_identical(p[names(shared)], shared)

  severity <- lancelot_parameters_severity(0.25, NULL)
  expect_identical(p[names(severity)], severity)

  observation <- lancelot_parameters_observation(1e6)
  expect_identical(p[names(observation)], observation)

  sens_and_spec <- lancelot_parameters_sens_and_spec()
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
      "strain_rel_p_ICU_D", "strain_rel_p_H_D",
      "strain_rel_p_W_D", "strain_rel_p_G_D",
      "strain_rel_p_icu",
      "strain_rel_p_hosp_if_sympt", "strain_rel_p_sympt", "N_tot_under15",
      "N_tot_15_24", "N_tot_25_49", "N_tot_50_64", "N_tot_65_79",
      "N_tot_80_plus", "N_tot_15_64", "N_tot_all", "N_tot_over25",
      "N_tot_react", "N_5_24_react", "N_25_34_react", "N_35_44_react",
      "N_45_54_react", "N_55_64_react", "N_65_plus_react", "I_A_transmission",
      "I_P_transmission", "I_C_1_transmission", "I_C_2_transmission",
      "n_groups", "initial_seed_size", "cross_immunity", "vacc_skip_from",
      "vacc_skip_to", "vacc_skip_dose", "vacc_skip_progression_rate_base",
      "vacc_skip_weight"))

  expect_equal(p$carehome_beds, sircovid_carehome_beds("uk"))
  expect_equal(p$carehome_residents, round(p$carehome_beds * 0.742))
  expect_equal(p$carehome_workers, p$carehome_residents)

  ## Can be slightly off due to rounding error
  expect_true(abs(sum(p$N_tot) - sum(p$population)) < 3)
  expect_length(p$N_tot, 19)
})


test_that("can compute severity for lancelot model", {
  severity <- lancelot_parameters_severity(0.25, NULL)

  expect_equal(
    severity$p_G_D_step, array(0.05, c(1, 19)))
  expect_equal(
    severity$p_star_step, array(0.2, c(1, 19)))
})


test_that("can input severity data for lancelot model", {
  data <- sircovid_parameters_severity(NULL)

  data$p_G_D[] <- 0
  data$p_star[] <- 0.5

  severity <- lancelot_parameters_severity(0.25, data)

  expect_equal(
    severity$p_G_D_step, array(0, c(1, 19)))
  expect_equal(
    severity$p_star_step, array(0.5, c(1, 19)))
})


test_that("can compute time-varying severity parameters for lancelot model", {
  dt <- 0.25

  p_G_D_date <- sircovid_date(c("2020-02-01", "2020-05-01"))
  p_G_D_value <- c(0.05, 0.1)
  p_G_D_CHR_value <- 0.4

  p_H_value <- 0.6
  p_H_CHR_date <- sircovid_date(c("2020-03-01", "2020-04-01"))
  p_H_CHR_value <- c(0.7, 0.6)

  severity <-
    lancelot_parameters_severity(dt, NULL,
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
    lancelot_parameters_severity(dt,
                                 p_C = list(date = 1,
                                            value = 0.3)),
    "As 'p_C' has a single 'value', expected NULL or missing 'date'")

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_ICU = list(date = c(1, 4, 5),
                                              value = c(0.2, 0.3))),
    "'date' and 'value' for 'p_ICU' must have the same length")

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_ICU_D = list(date = c(1, 4),
                                                value = c(-1, 0.3))),
    "'p_ICU_D' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_W_D = list(date = c(1, 4),
                                              value = c(0.2, 3))),
    "'p_W_D' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_H_CHR = list(date = 1,
                                                value = 0.3)),
    "As 'p_H_CHR' has a single 'value', expected NULL or missing 'date'")

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_G_D_CHR = list(date = c(1, 4, 5),
                                                  value = c(0.2, 0.3))),
    "'date' and 'value' for 'p_G_D_CHR' must have the same length")

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_H_CHR = list(date = c(1, 4),
                                                value = c(-1, 0.3))),
    "'p_H_CHR' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    lancelot_parameters_severity(dt,
                                 p_G_D_CHR = list(date = c(1, 4),
                                                  value = c(0.2, 3))),
    "'p_G_D_CHR' must lie in [0, 1]", fixed = TRUE)

})


test_that("can compute time-varying progression parameters for lancelot
          model", {
            dt <- 0.25

            gamma_H_R_value <- 0.3
            gamma_H_D_date <- sircovid_date(c("2020-02-01", "2020-05-01"))
            gamma_H_D_value <- c(0.2, 0.5)

            progression <-
              lancelot_parameters_progression(dt,
                                              gamma_H_D =
                                                list(date = gamma_H_D_date,
                                                     value = gamma_H_D_value),
                                              gamma_H_R =
                                                list(value = gamma_H_R_value)
              )

            gamma_H_D_step <-
              sircovid_parameters_piecewise_linear(gamma_H_D_date,
                                                   gamma_H_D_value, dt)
            expect_equal(progression$gamma_H_D_step, gamma_H_D_step)
            expect_equal(progression$gamma_H_R_step, gamma_H_R_value)

            expect_error(
              lancelot_parameters_progression(dt,
                                              gamma_E = list(date = 1,
                                                             value = 3)),
              "'gamma_E' has a single 'value', expected NULL or missing 'date'")

            expect_error(
              lancelot_parameters_progression(dt,
                                              gamma_ICU_pre =
                                                list(date = c(1, 4, 5),
                                                     value = c(2, 3))),
              "'date' and 'value' for 'gamma_ICU_pre' must have the same length"
            )

            expect_error(
              lancelot_parameters_progression(dt,
                                              gamma_H_D =
                                                list(date = c(1, 4),
                                                     value = c(-2, 3))),
              "'gamma_H_D' must have only non-negative values")


          })


test_that("Can compute transmission matrix for lancelot model", {
  m <- lancelot_transmission_matrix(0.1, 4e-5, 5e-4, "uk")
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
  p1 <- lancelot_parameters_observation(1e6)
  p2 <- lancelot_parameters_observation(1e4)
  expect_setequal(names(p1), names(p2))
  v <- setdiff(names(p1), "exp_noise")
  expect_mapequal(p1[v], p2[v])
  expect_equal(p1$exp_noise, 1e6)
  expect_equal(p2$exp_noise, 1e4)
})


test_that("lancelot_index identifies ICU and D_tot in real model", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  mod <- lancelot$new(p, 0, 10)
  info <- mod$info()
  index <- lancelot_index(info)

  expect_equal(
    names(index$run),
    c("time", "icu", "general", "deaths_carehomes_inc", "deaths_comm_inc",
      "deaths_comm_0_49_inc", "deaths_comm_50_54_inc", "deaths_comm_55_59_inc",
      "deaths_comm_60_64_inc", "deaths_comm_65_69_inc", "deaths_comm_70_74_inc",
      "deaths_comm_75_79_inc", "deaths_comm_80_plus_inc",
      "deaths_hosp_inc", "deaths_hosp_0_49_inc", "deaths_hosp_50_54_inc",
      "deaths_hosp_55_59_inc", "deaths_hosp_60_64_inc", "deaths_hosp_65_69_inc",
      "deaths_hosp_70_74_inc", "deaths_hosp_75_79_inc",
      "deaths_hosp_80_plus_inc", "admitted_inc", "all_admission_0_9_inc",
      "all_admission_10_19_inc", "all_admission_20_29_inc",
      "all_admission_30_39_inc", "all_admission_40_49_inc",
      "all_admission_50_59_inc", "all_admission_60_69_inc",
      "all_admission_70_79_inc", "all_admission_80_plus_inc", "diagnoses_inc",
      "sero_pos_1", "sero_pos_2", "sympt_cases_inc",
      "sympt_cases_non_variant_inc", "sympt_cases_over25_inc",
      "sympt_cases_under15_inc", "sympt_cases_15_24_inc",
      "sympt_cases_25_49_inc", "sympt_cases_50_64_inc",
      "sympt_cases_65_79_inc", "sympt_cases_80_plus_inc",
      "sympt_cases_non_variant_over25_inc", "react_pos",
      "react_5_24_pos", "react_25_34_pos", "react_35_44_pos",
      "react_45_54_pos", "react_55_64_pos", "react_65_plus_pos"))

  expect_equal(index$run[["time"]],
               which(names(info$index) == "time"))
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
  expect_equal(lancelot_index_workers(), 6:13) })


test_that("lancelot_index is properly named", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  mod <- lancelot$new(p, 0, 10)
  info <- mod$info()
  index <- lancelot_index(info)
  expect_false(any(is.na(names(index))))
})


test_that("Can compute initial conditions", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 10)
  mod <- lancelot$new(p, 0, 10)
  info <- mod$info()

  initial <- lancelot_initial(info, 10, p)
  initial_y <- mod$transform_variables(initial)

  expect_equal(initial_y$N_tot_sero_1, sum(p$N_tot))
  expect_equal(initial_y$N_tot_sero_2, sum(p$N_tot))
  expect_equal(initial_y$N_tot_PCR, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S),
               p$N_tot)

  ## 48 here, derived from;
  ## * 38 (S + N_tot)
  ## * 1 (react_pos)
  ## * 3 (N_tot_sero_1 + N_tot_sero_2 + N_tot_PCR)
  expect_equal(sum(initial != 0), 43)
})


test_that("Can control the seeding", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "london",
                           initial_seed_size = 50)
  expect_equal(p$initial_seed_size, 50)

  mod <- lancelot$new(p, 0, 10)
  info <- mod$info()

  initial <- lancelot_initial(info, 10, p)

  initial_y <- mod$transform_variables(initial)

  expect_equal(initial_y$N_tot_sero_1, sum(p$N_tot))
  expect_equal(initial_y$N_tot_sero_2, sum(p$N_tot))
  expect_equal(initial_y$N_tot_PCR, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S) + drop(initial_y$I_A),
               p$N_tot)

  ## 48 here, derived from;
  ## * 38 (S + N_tot)
  ## * 1 (prob_strain)
  ## * 3 (N_tot_sero_1 + N_tot_sero_2 + N_tot_PCR)
  expect_equal(sum(initial != 0), 43)
})

## TODO: Ed - you had said that you had ideas for some more systematic
## testing here.  This will also get easier to do if we move to having
## a function generator given a data set.
test_that("lancelot_compare combines likelihood correctly", {
  state <- rbind(
    time = 45,
    icu = 10:15,
    general = 20:25,
    deaths_carehomes_inc = 2:7,
    deaths_comm_inc = 1:6,
    deaths_comm_0_49_inc = 1:6,
    deaths_comm_50_54_inc = 1:6,
    deaths_comm_55_59_inc = 2:7,
    deaths_comm_60_64_inc = 2:7,
    deaths_comm_65_69_inc = 3:8,
    deaths_comm_70_74_inc = 3:8,
    deaths_comm_75_79_inc = 4:9,
    deaths_comm_80_plus_inc = 4:9,
    deaths_hosp_inc = 3:8,
    deaths_hosp_0_49_inc = 1:6,
    deaths_hosp_50_54_inc = 1:6,
    deaths_hosp_55_59_inc = 2:7,
    deaths_hosp_60_64_inc = 2:7,
    deaths_hosp_65_69_inc = 3:8,
    deaths_hosp_70_74_inc = 3:8,
    deaths_hosp_75_79_inc = 4:9,
    deaths_hosp_80_plus_inc = 4:9,
    admitted_inc = 50:55,
    diagnoses_inc = 60:65,
    all_admission_0_9_inc = 1:6,
    all_admission_10_19_inc = 1:6,
    all_admission_20_29_inc = 1:6,
    all_admission_30_39_inc = 1:6,
    all_admission_40_49_inc = 9:14,
    all_admission_50_59_inc = 12:17,
    all_admission_60_69_inc = 17:22,
    all_admission_70_79_inc = 20:25,
    all_admission_80_plus_inc = 27:32,
    sero_pos_1 = 4:9,
    sero_pos_2 = 14:19,
    sympt_cases_inc = 100:105,
    sympt_cases_non_variant_inc = 70:75,
    sympt_cases_over25_inc = 80:85,
    sympt_cases_under15_inc = 5:10,
    sympt_cases_15_24_inc = 5:10,
    sympt_cases_25_49_inc = 19:24,
    sympt_cases_50_64_inc = 19:24,
    sympt_cases_65_79_inc = 19:24,
    sympt_cases_80_plus_inc = 19:24,
    sympt_cases_non_variant_over25_inc = 60:65,
    react_pos = 2:7,
    react_5_24_pos = 1:6,
    react_25_34_pos = 1:6,
    react_35_44_pos = 1:6,
    react_45_54_pos = 1:6,
    react_55_64_pos = 1:6,
    react_65_plus_pos = 1:6)
  observed <- list(
    icu = 13,
    general = 23,
    hosp = 36,
    deaths_carehomes = 4,
    deaths_hosp = 5,
    deaths_hosp_0_49 = 1,
    deaths_hosp_50_54 = 2,
    deaths_hosp_55_59 = 3,
    deaths_hosp_60_64 = 4,
    deaths_hosp_65_69 = 5,
    deaths_hosp_70_74 = 6,
    deaths_hosp_75_79 = 7,
    deaths_hosp_80_plus = 8,
    deaths_comm = 3,
    deaths_comm_0_49 = 1,
    deaths_comm_50_54 = 2,
    deaths_comm_55_59 = 3,
    deaths_comm_60_64 = 4,
    deaths_comm_65_69 = 5,
    deaths_comm_70_74 = 6,
    deaths_comm_75_79 = 7,
    deaths_comm_80_plus = 8,
    deaths = 8,
    deaths_non_hosp = 6,
    admitted = 53,
    diagnoses = 63,
    all_admission = 116,
    all_admission_0_9 = 2,
    all_admission_10_19 = 4,
    all_admission_20_29 = 5,
    all_admission_30_39 = 5,
    all_admission_40_49 = 10,
    all_admission_50_59 = 15,
    all_admission_60_69 = 20,
    all_admission_70_79 = 25,
    all_admission_80_plus = 30,
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
    pillar2_under15_cases = 8,
    pillar2_15_24_cases = 8,
    pillar2_25_49_cases = 20,
    pillar2_50_64_cases = 20,
    pillar2_65_79_cases = 20,
    pillar2_80_plus_cases = 20,
    pillar2_under15_pos = 8,
    pillar2_15_24_pos = 8,
    pillar2_25_49_pos = 20,
    pillar2_50_64_pos = 20,
    pillar2_65_79_pos = 20,
    pillar2_80_plus_pos = 20,
    pillar2_under15_tot = 160,
    pillar2_15_24_tot = 160,
    pillar2_25_49_tot = 400,
    pillar2_50_64_tot = 400,
    pillar2_65_79_tot = 400,
    pillar2_80_plus_tot = 400,
    react_pos = 3,
    react_tot = 500,
    react_5_24_pos = 1,
    react_5_24_tot = 50,
    react_25_34_pos = 1,
    react_25_34_tot = 50,
    react_35_44_pos = 1,
    react_35_44_tot = 100,
    react_45_54_pos = 1,
    react_45_54_tot = 100,
    react_55_64_pos = 1,
    react_55_64_tot = 100,
    react_65_plus_pos = 2,
    react_65_plus_tot = 100,
    strain_non_variant = 40,
    strain_tot = 50,
    strain_over25_non_variant = 20,
    strain_over25_tot = 25)
  date <- sircovid_date("2020-01-01")
  pars <- lancelot_parameters(date, "uk", exp_noise = Inf)

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
  nms_pillar2_under15 <- c("pillar2_under15_pos", "pillar2_under15_tot")
  nms_pillar2_15_24 <- c("pillar2_15_24_pos", "pillar2_15_24_tot")
  nms_pillar2_25_49 <- c("pillar2_25_49_pos", "pillar2_25_49_tot")
  nms_pillar2_50_64 <- c("pillar2_50_64_pos", "pillar2_50_64_tot")
  nms_pillar2_65_79 <- c("pillar2_65_79_pos", "pillar2_65_79_tot")
  nms_pillar2_80_plus <- c("pillar2_80_plus_pos", "pillar2_80_plus_tot")
  nms_react <- c("react_pos", "react_tot")
  nms_react_5_24 <- c("react_5_24_pos", "react_5_24_tot")
  nms_react_25_34 <- c("react_25_34_pos", "react_25_34_tot")
  nms_react_35_44 <- c("react_35_44_pos", "react_35_44_tot")
  nms_react_45_54 <- c("react_45_54_pos", "react_45_54_tot")
  nms_react_55_64 <- c("react_55_64_pos", "react_55_64_tot")
  nms_react_65_plus <- c("react_65_plus_pos", "react_65_plus_tot")
  nms_strain <- c("strain_non_variant", "strain_tot")
  nms_strain_over25 <- c("strain_over25_non_variant", "strain_over25_tot")
  parts <- c(as.list(setdiff(names(observed),
                             c(nms_sero_1, nms_sero_2, nms_pillar2,
                               nms_pillar2_over25, nms_pillar2_under15,
                               nms_pillar2_15_24, nms_pillar2_25_49,
                               nms_pillar2_50_64, nms_pillar2_65_79,
                               nms_pillar2_80_plus, nms_react, nms_react_5_24,
                               nms_react_25_34, nms_react_35_44,
                               nms_react_45_54, nms_react_55_64,
                               nms_react_65_plus, nms_strain,
                               nms_strain_over25))),
             list(nms_sero_1), list(nms_sero_2), list(nms_pillar2),
             list(nms_pillar2_over25), list(nms_pillar2_under15),
             list(nms_pillar2_15_24), list(nms_pillar2_25_49),
             list(nms_pillar2_50_64), list(nms_pillar2_65_79),
             list(nms_pillar2_80_plus), list(nms_react),
             list(nms_react_5_24), list(nms_react_25_34),
             list(nms_react_35_44), list(nms_react_45_54),
             list(nms_react_55_64), list(nms_react_65_plus),
             list(nms_strain), list(nms_strain_over25))

  ll_parts <- lapply(parts, function(x)
    lancelot_compare(state, observed_keep(x), pars))

  ## Extremely light testing, though this has already flushed out some
  ## issues
  expect_vector_equal(lengths(ll_parts), 6)
  expect_equal(
    lancelot_compare(state, observed, pars),
    rowSums(do.call(cbind, ll_parts)))

  ## Test that there are non-zero values for each log-likelihood part.
  ## This helps make sure all parts contribute to the log-likelihood.
  expect_true(all(sapply(ll_parts, function(x) any(x != 0))))


  ## Check that weekend effect parameters work as expected. For each day
  ## use same state values except time
  state <- state[, 1, drop = FALSE]
  time <- seq(20, 26, 1)
  age_nms <- c("_under15", "_15_24", "_25_49",
               "_50_64", "_65_79", "_80_plus")

  for (i in age_nms) {
    pars[[paste0("p_NC_weekend", i)]] <- pars[[paste0("p_NC", i)]]
    pars[[paste0("phi_pillar2_cases_weekend", i)]] <-
      pars[[paste0("phi_pillar2_cases", i)]]
  }

  helper <- function(parameter, age_nm) {

    pars[[paste0(parameter, "_weekend", age_nm)]] <-
      0.5 * pars[[paste0(parameter, age_nm)]]

    ll_time <- vnapply(time, function(x) {
      state["time", ] <- x
      lancelot_compare(state, observed, pars)
    })

    ## First 5 days are weekdays, last 2 are weekend
    expect_equal(grepl("^S", weekdays(sircovid_date_as_date(time))),
                 c(rep(FALSE, 5), rep(TRUE, 2)))
    ## Weekdays should yield the same values. Weekends should yield the same
    ## values, but different to weekdays
    expect_equal(ll_time[2:5], rep(ll_time[1], 4))
    expect_equal(ll_time[6], ll_time[7])
    expect_true(ll_time[1] != ll_time[6])
  }

  lapply(age_nms, function(x) helper("p_NC", x))
  lapply(age_nms, function(x) helper("phi_pillar2_cases", x))

})


test_that("lancelot_population prevents negative populations", {
  population <- rep(100, 17)
  expect_error(
    lancelot_population(population, 1000, 0),
    "Not enough population to be care workers")
  expect_error(
    lancelot_population(population, 0, 1000),
    "Not enough population to meet care home occupancy")
})


test_that("lancelot_population preserves population", {
  population <- c(949, 989, 1030, 995, 993, 985, 965, 1082, 1042, 960,
                  980, 1004, 934, 1049, 971, 1020, 937)
  res <- lancelot_population(population, 120, 200)
  expect_equal(sum(res), sum(population))
  expect_equal(res[18:19], c(120, 200))
  expect_equal(res[1:5], population[1:5])
  expect_vector_lt(res[6:17], population[6:17])
})


test_that("lancelot_index switches work as expected", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  mod <- lancelot$new(p, 0, 5, seed = 1L)
  info <- mod$info()

  index <- lancelot_index(info)
  expect_equal(
    unname(index$state[!names(index$state) %in% names(index$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))

  index_rt <- lancelot_index(info, rt = FALSE)
  expect_equal(
    unname(index_rt$state[!names(index_rt$state) %in% names(index_rt$run)]),
     c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))


  index_cum_admit <- lancelot_index(info, cum_admit = FALSE)
  expect_equal(
    unname(index_cum_admit$state[!names(index_cum_admit$state) %in%
    names(index_cum_admit$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))

  index_diagnoses_admitted <- lancelot_index(info, diagnoses_admitted = FALSE)
  expect_equal(
    unname(
      index_diagnoses_admitted$state[!names(index_diagnoses_admitted$state)
      %in% names(index_diagnoses_admitted$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))

  index_cum_infections_disag <- lancelot_index(info,
   cum_infections_disag = FALSE)
  expect_equal(
    unname(
      index_cum_infections_disag$state[!names(index_cum_infections_disag$state)
      %in% names(index_cum_infections_disag$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))

  index_cum_n_vaccinated <- lancelot_index(info, cum_n_vaccinated = FALSE)
  expect_equal(
    unname(
      index_cum_n_vaccinated$state[!names(index_cum_n_vaccinated$state)
      %in% names(index_cum_n_vaccinated$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$D,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))

  index_inf_per_strain <-
    lancelot_index(info, infections_inc_per_strain = FALSE)
  expect_equal(
    unname(
      index_inf_per_strain$state[!names(index_inf_per_strain$state)
                                   %in% names(index_inf_per_strain$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$D_hosp))

  index_D <- lancelot_index(info, D_all = FALSE)
  expect_equal(
    unname(
      index_D$state[!names(index_D$state)
      %in% names(index_D$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D_hosp,
      info$index$infections_inc_per_strain))

  index_D_hosp <- lancelot_index(info, D_hosp = FALSE)
  expect_equal(
    unname(
      index_D_hosp$state[!names(index_D_hosp$state)
                         %in% names(index_D_hosp$run)]),
    c(info$index$hosp_tot,
      info$index$cum_admit_conf,
      info$index$cum_new_conf,
      info$index$D_tot,
      info$index$D_inc,
      info$index$cum_infections,
      info$index$infections_inc,
      info$index$effective_susceptible,
      info$index$S,
      info$index$R,
      info$index$prob_strain,
      info$index$cum_admit_by_age,
      info$index$diagnoses_admitted,
      info$index$cum_infections_disag,
      info$index$cum_n_vaccinated,
      info$index$D,
      info$index$infections_inc_per_strain))
})


test_that("lancelot_check_data requires consistent deaths", {
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))
  data$deaths_hosp <- data$deaths
  expect_error(
    lancelot_check_data(data),
    paste("Deaths are not consistently split into total vs",
          "hospital/non-hospital or hospital/care homes/community"))

  data$deaths_non_hosp <- data$deaths
  data$deaths_carehomes <- data$deaths
  data$deaths <- NA_real_
  expect_error(
    lancelot_check_data(data),
    paste("Non-hospital deaths are not consistently split into total vs",
          "care homes/community"))
})


test_that("lancelot_check_data disallows double fitting to hospital deaths", {
  ## Specifically, disallow double fitting to aggregated and
  ## disaggregated hospital deaths
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))

  data$deaths <- NA
  data$deaths_hosp <- 10

  expect_silent(lancelot_check_data(data))

  ## Allow age disaggregated and aggregated on different days
  data$deaths_hosp[11:20] <- NA
  data$deaths_hosp_60_64[11:20] <- 3
  data$deaths_hosp_80_plus[11:20] <- 7

  expect_silent(lancelot_check_data(data))

  ## Throw error
  data$deaths_hosp_0_49 <- 1
  data$deaths_hosp_65_69 <- 4
  expect_error(
    lancelot_check_data(data),
    paste0("Cannot fit to all ages aggregated for hospital deaths if fitting ",
           "to any sub-groups"))

  ## Add populations and throw regional error
  data2 <- cbind(
    rbind(data, data),
    region = factor(rep(c("london", "south_east"), each = nrow(data))))

  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to all ages aggregated for hospital deaths if")
})


test_that("lancelot_check_data disallows double fitting to community deaths", {
  ## Specifically, disallow double fitting to aggregated and
  ## disaggregated community deaths
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))

  data$deaths <- NA
  data$deaths_comm <- 10

  expect_silent(lancelot_check_data(data))

  ## Allow age disaggregated and aggregated on different days
  data$deaths_comm[11:20] <- NA
  data$deaths_comm_60_64[11:20] <- 3
  data$deaths_comm_80_plus[11:20] <- 7

  expect_silent(lancelot_check_data(data))

  ## Throw error
  data$deaths_comm_0_49 <- 1
  data$deaths_comm_65_69 <- 4
  expect_error(
    lancelot_check_data(data),
    paste0("Cannot fit to all ages aggregated for community deaths if fitting ",
           "to any sub-groups"))

  ## Add populations and throw regional error
  data2 <- cbind(
    rbind(data, data),
    region = factor(rep(c("london", "south_east"), each = nrow(data))))

  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to all ages aggregated for community deaths if")
})


test_that("lancelot_check_data disallows double fitting to admissions", {
  ## Specifically, disallow double fitting to aggregated and
  ## age-specific hospital admissions
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))

  expect_silent(lancelot_check_data(data))

  ## Expect no error
  n <- nrow(data)
  data[1:n - 1, "all_admission_60_69"] <- 1
  data[n, "all_admission"] <- 4
  expect_silent(lancelot_check_data(data))

  ## Throw error
  data$all_admission_60_69 <- 1
  data$all_admission <- 4
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to admissions by age and aggregate together!")

  ## Add populations and throw regional error
  data2 <- cbind(
    rbind(data, data),
    region = factor(rep(c("london", "south_east"), each = nrow(data))))
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to admissions by age and aggregate together!")
})


test_that("lancelot_check_data disallows double fitting to REACT", {
  ## Specifically, disallow double fitting to aggregated and
  ## age-specific hospital admissions
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))

  expect_silent(lancelot_check_data(data))

  ## Expect no error
  n <- nrow(data)
  data[1:n - 1, "react_25_34_pos"] <- 1
  data[n, "react_pos"] <- 4
  expect_silent(lancelot_check_data(data))

  ## Throw error
  data$react_25_34_pos <- 1
  data$react_pos <- 4
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to REACT by age and aggregate together!")

  ## Add populations and throw regional error
  data2 <- cbind(
    rbind(data, data),
    region = factor(rep(c("london", "south_east"), each = nrow(data))))
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to REACT by age and aggregate together!")
})


test_that("lancelot_check_data prevents pillar 2 double fitting", {
  ## does not allow more than one pillar 2, strain data stream or
  ## pillar 2 double fitting
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))

  ## Example that should not cause error
  data$pillar2_pos <- 3
  data$pillar2_tot <- 6
  data$strain_non_variant <- 8
  data$strain_tot <- 10
  expect_silent(lancelot_check_data(data))

  data$pillar2_pos <- 3
  data$pillar2_tot <- 6
  data$pillar2_cases <- 5
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to pillar 2 cases and positivity together")

  data$pillar2_pos <- 3
  data$pillar2_tot <- 6
  data$pillar2_cases <- NA
  data$pillar2_25_49_cases <- 5
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to pillar 2 cases and positivity together")

  data$pillar2_pos <- 3
  data$pillar2_tot <- 6
  data$pillar2_25_49_cases <- NA
  data$pillar2_under15_pos <- 3
  data$pillar2_under15_tot <- 6
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to all ages aggregated for pillar 2 if fitting to")

  data$pillar2_pos <- NA
  data$pillar2_tot <- NA
  data$pillar2_over25_cases <- 5
  data$pillar2_65_79_cases <- 4
  data$pillar2_under15_pos <- NA
  data$pillar2_under15_tot <- NA
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to over 25s for pillar 2 if fitting to any over 25")

  data$pillar2_over25_cases <- NA
  data$pillar2_65_79_cases <- NA
  data$strain_non_variant <- 8
  data$strain_tot <- 10
  data$strain_over25_non_variant <- 6
  data$strain_over25_tot <- 9
  expect_error(
    lancelot_check_data(data),
    "Cannot fit to more than one strain data stream")
})


test_that("the lancelot sircovid model has 19 groups", {
  expect_equal(lancelot_n_groups(), 19)
})


test_that("lancelot_check_data prevents pillar 2 double fitting (nested)", {
  ## does not allow more than one pillar 2, strain data stream or
  ## pillar 2 double fitting (for grouped data)
  data <- lancelot_simple_data(read_csv(sircovid_file("extdata/example.csv")))
  data2 <- cbind(
    rbind(data, data),
    region = factor(rep(c("london", "south_east"), each = nrow(data))))

  ## Example that should not cause error
  data2$pillar2_pos <- 3
  data2$pillar2_tot <- 6
  expect_silent(lancelot_check_data(data2))

  ## Error if one region has multiple streams
  data2$pillar2_pos <- 3
  data2$pillar2_tot <- 6
  data2$pillar2_cases <- 5
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to pillar 2 cases and positivity together")

  data2$pillar2_pos <- 3
  data2$pillar2_tot <- 6
  data2$pillar2_cases <- NA
  data2$pillar2_25_49_cases <- 5
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to pillar 2 cases and positivity together")

  data2$pillar2_pos <- 3
  data2$pillar2_tot <- 6
  data2$pillar2_25_49_cases <- NA
  data2$pillar2_under15_pos <- 3
  data2$pillar2_under15_tot <- 6
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to all ages aggregated for pillar 2 if fitting")

  data2$pillar2_pos <- NA
  data2$pillar2_tot <- NA
  data2$pillar2_over25_cases <- 5
  data2$pillar2_65_79_cases <- 4
  data2$pillar2_under15_pos <- NA
  data2$pillar2_under15_tot <- NA
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to over 25s for pillar 2 if fitting to any over")

  data2$pillar2_over25_cases <- NA
  data2$pillar2_65_79_cases <- NA
  data2$strain_non_variant <- 8
  data2$strain_tot <- 10
  data2$strain_over25_non_variant <- 6
  data2$strain_over25_tot <- 9
  expect_error(
    lancelot_check_data(data2),
    "london: Cannot fit to more than one strain data stream")

  ## Don't error if two regions have different streams
  i1 <- seq_len(nrow(data))
  i2 <- i1 + nrow(data)
  data2[i1, "pillar2_cases"] <- 3
  data2[i2, "pillar2_pos"] <- 3
  data2[i2, "pillar2_tot"] <- 5
  data2[i1, "strain_non_variant"] <- NA
  data2[i1, "strain_tot"] <- NA
  data2[i2, "strain_over25_non_variant"] <- NA
  data2[i2, "strain_over25_tot"] <- NA
  expect_silent(lancelot_check_data(data2))
})


test_that("lancelot check severity works as expected", {
  ## errors if params missing
  expect_error(
    lancelot_check_severity(list(n_groups = 19)),
    "Parameter 'rel_p_sympt' is missing"
  )
  expect_error(
    lancelot_check_severity(list(n_groups = 19, rel_p_sympt = 1)),
    "Parameter 'p_C_step' is missing"
  )

  ## errors if rel_ is not 1 or 4 cols
  expect_error(lancelot_check_severity(list(
    n_groups = 19,
    rel_p_sympt = array(1, c(19, 3, 3)),
    p_C_step = matrix(1, 1, 19)
  )), "1 or 4 columns")

  ## no error if 1 col
  expect_error(lancelot_check_severity(list(
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
  expect_equal(lancelot_check_severity(p), p)

  ## check errors when expected - we only check upper bound as it's safe to
  ##  assume negative probs won't be given in practice (but uses the same )

  ## rel_p_sympt Inf
  for (i in seq_along(steps)) {
    p[[rels[[i]]]][] <- -1
    expect_error(
      lancelot_check_severity(p),
      sprintf("%s * %s is not in [0, 1] for group 1", rels[[i]], steps[[i]]),
      fixed = TRUE
    )
    p[[rels[[i]]]][] <- 1
  }
})


test_that("Can input population data", {

  pop <- sircovid_population("uk")
  p <- lancelot_parameters(1, "Unknown", population = pop, carehome_beds = 0)
  expect_equal(p$population, pop)
  expect_equal(p$N_tot, c(pop, 0, 0))

  expect_error(
    lancelot_parameters(1, "Unknown", population = pop[-1L],
                        carehome_beds = 0),
    "If population is specified it must be a vector of length 17")

  pop[1L] <- 0.5
  expect_error(
    lancelot_parameters(1, "Unknown", population = pop, carehome_beds = 0),
    "'population' must be an integer")

  pop[1L] <- -1
  expect_error(
    lancelot_parameters(1, "Unknown", population = pop, carehome_beds = 0),
    "'population' must have only non-negative values")
})

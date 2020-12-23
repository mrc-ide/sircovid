context("carehomes (support)")

test_that("carehomes progression parameters", {
  p <- carehomes_parameters_progression()
  expect_setequal(
    names(p),
    c("s_E", "s_asympt", "s_sympt", "s_comm_D", "s_hosp_D", "s_hosp_R",
      "s_ICU_D", "s_ICU_S_R", "s_ICU_S_D", "s_triage", "s_stepdown_D",
      "s_stepdown_R", "s_R_pos", "s_PCR_pos", "s_PCR_pre", "gamma_E",
      "gamma_asympt", "gamma_sympt", "gamma_comm_D", "gamma_hosp_D",
      "gamma_hosp_R", "gamma_ICU_D", "gamma_ICU_S_R", "gamma_ICU_S_D",
      "gamma_triage", "gamma_stepdown_D", "gamma_stepdown_R", "gamma_R_pos",
      "gamma_R_pre_1", "gamma_R_pre_2", "gamma_test", "gamma_PCR_pos",
      "gamma_PCR_pre"))

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
    c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt",
      "vaccine_progression_rate_base", "vaccine_population_reluctant",
      "vaccine_daily_doses"))
  expect_equal(nrow(p$rel_susceptibility), n_groups)
  expect_equal(ncol(p$rel_susceptibility), 1)
  expect_equal(nrow(p$vaccine_progression_rate_base), n_groups)
  expect_equal(ncol(p$vaccine_progression_rate_base), 1)

  # test when more vaccinated categories than default
  rel_susceptibility <- c(1, 0.2, 0.1, 0.4)
  rel_p_sympt <- c(1, 0.75, 0.5, 0.75)
  rel_p_hosp_if_sympt <- c(1, 0.8, 0.6, 0.9)
  vaccine_progression_rate <- c(0, 1, 1, 1)
  p <- carehomes_parameters_vaccination(ntot,
                                        rel_susceptibility = rel_susceptibility,
                                        rel_p_sympt = rel_p_sympt,
                                        rel_p_hosp_if_sympt =
                                          rel_p_hosp_if_sympt,
                                        vaccine_progression_rate =
                                          vaccine_progression_rate)
  expect_setequal(
    names(p),
    c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt",
      "vaccine_progression_rate_base", "vaccine_population_reluctant",
      "vaccine_daily_doses"))
  expect_equal(nrow(p$rel_susceptibility), n_groups)
  expect_equal(ncol(p$rel_susceptibility), length(rel_susceptibility))
  expect_equal(nrow(p$rel_p_sympt), n_groups)
  expect_equal(ncol(p$rel_p_sympt), length(rel_p_sympt))
  expect_equal(nrow(p$rel_p_hosp_if_sympt), n_groups)
  expect_equal(ncol(p$rel_p_hosp_if_sympt), length(rel_p_hosp_if_sympt))
  expect_equal(nrow(p$vaccine_progression_rate_base), n_groups)
  expect_equal(ncol(p$vaccine_progression_rate_base),
               length(vaccine_progression_rate))
  msg1 <- "rel_susceptibility, rel_p_sympt, rel_p_hosp_if_sympt"
  msg2 <- "should have the same dimension"
  expect_error(
    carehomes_parameters_vaccination(ntot,
                                     rel_susceptibility = 1,
                                     rel_p_sympt = c(1, 0.5, 0.25),
                                     rel_p_hosp_if_sympt = c(1, 0.1)),
    paste(msg1, msg2))
  expect_error(carehomes_parameters_vaccination(ntot,
                                                rel_susceptibility = c(1, 1),
                                                rel_p_sympt = c(1, 0.5, 0.25),
                                                rel_p_hosp_if_sympt = 1),
               paste(msg1, msg2))
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

  progression <- carehomes_parameters_progression()
  expect_identical(p[names(progression)], progression)

  vaccination <- carehomes_parameters_vaccination(
    p$N_tot, p$rel_susceptibility, p$rel_p_sympt, p$rel_p_hosp_if_sympt,
    p$vaccine_progression_rate_base)
  expect_identical(p[names(vaccination)], vaccination)
  
  strain <- carehomes_parameters_strain(p$strain_transmission,
                                        strain_seed_date = NULL,
                                        strain_seed_value = NULL, 
                                        dt = 1/4)

  waning <- carehomes_parameters_waning(0)
  expect_identical(p[names(waning)], waning)

  shared <- sircovid_parameters_shared(date, "uk", NULL, NULL)
  expect_identical(p[names(shared)], shared)

  severity <- carehomes_parameters_severity(NULL, 0.7)
  expect_identical(p[names(severity)], severity)

  expect_equal(
    p$observation,
    carehomes_parameters_observation(1e6))
  expect_equal(p$N_tot_15_64, sum(p$N_tot[4:13]))

  extra <- setdiff(names(p),
                   c("m", "observation",
                     names(shared), names(progression), names(severity),
                     names(strain), names(vaccination), names(waning),
                     "model_pcr_and_serology_user"))
  expect_setequal(
    extra,
    c("N_tot", "carehome_beds", "carehome_residents", "carehome_workers",
      "sero_specificity", "sero_sensitivity", "N_tot_15_64",
      "pillar2_specificity", "pillar2_sensitivity", "react_specificity",
      "react_sensitivity", "prop_noncovid_sympt", "psi_death_ICU",
      "p_death_ICU_step", "psi_death_hosp_D", "p_death_hosp_D_step",
      "psi_death_stepdown", "p_death_stepdown_step", "psi_hosp_sympt",
      "p_hosp_sympt_step", "psi_death_comm", "p_death_comm_step",
      "psi_ICU_hosp", "p_ICU_hosp_step", "psi_admit_conf", "p_admit_conf_step",
      "n_groups"))

  expect_equal(p$carehome_beds, sircovid_carehome_beds("uk"))
  expect_equal(p$carehome_residents, round(p$carehome_beds * 0.742))
  expect_equal(p$carehome_workers, p$carehome_residents)

  ## Can be slightly off due to rounding error
  expect_true(abs(sum(p$N_tot) - sum(p$population)) < 3)
  expect_length(p$N_tot, 19)
})


test_that("can compute severity for carehomes model", {
  population <- sircovid_population("uk")
  severity <- carehomes_parameters_severity(NULL, 0.7)
  expect_true(all(lengths(severity) == 19))
  expect_setequal(names(severity), names(sircovid_parameters_severity(NULL)))

  expect_true(
    all(severity$p_serocoversion == severity$p_serocoversion[[1]]))
  expect_equal(
    severity$p_death_comm, rep(c(0, 0.7), c(18, 1)))
  expect_equal(
    severity$p_admit_conf, rep(0.2, 19))
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
    c("icu", "general", "deaths_comm", "deaths_hosp", "admitted", "new",
      "sero_pos", "sympt_cases", "sympt_cases_over25", "react_pos"))

  expect_equal(index$run[["icu"]],
               which(names(info$index) == "I_ICU_tot"))
  expect_equal(index$run[["general"]],
               which(names(info$index) == "general_tot"))
  expect_equal(index$run[["deaths_comm"]],
               which(names(info$index) == "D_comm_tot"))
  expect_equal(index$run[["deaths_hosp"]],
               which(names(info$index) == "D_hosp_tot"))
  expect_equal(index$run[["admitted"]],
               which(names(info$index) == "cum_admit_conf"))
  expect_equal(index$run[["new"]],
               which(names(info$index) == "cum_new_conf"))
  expect_equal(index$run[["sero_pos"]],
               which(names(info$index) == "sero_pos"))
  expect_equal(index$run[["sympt_cases"]],
               which(names(info$index) == "cum_sympt_cases"))
  expect_equal(index$run[["sympt_cases_over25"]],
               which(names(info$index) == "cum_sympt_cases_over25"))
  expect_equal(index$run[["react_pos"]],
               which(names(info$index) == "react_pos"))
})


test_that("carehome worker index is sensible", {
  expect_equal(carehomes_index_workers(), 6:13)
})


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

  expect_equal(initial_y$N_tot3, sum(p$N_tot))
  expect_equal(initial_y$N_tot2, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S) + drop(initial_y$I_asympt),
               p$N_tot)
  expect_equal(drop(initial_y$I_asympt),
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$R_pre[, 1, 1, ],
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$PCR_pos[, 1, 1, ],
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$react_pos, 10)

  ## 42 here, derived from;
  ## * 19 (S)
  ## * 19 (N_tot)
  ## * 4 values as N_tot2 + N_tot3 + I_asympt[4] + R_pre[4] + PCR_pos[4]
  expect_equal(sum(initial$state != 0), 44)
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
    deaths_comm = 1:6,
    deaths_hosp = 3:8,
    admitted = 50:55,
    new = 60:65,
    sero_pos = 4:9,
    sympt_cases = 100:105,
    sympt_cases_over25 = 80:85,
    react_pos = 2:7)
  prev_state <- array(1, dim(state), dimnames = dimnames(state))
  observed <- list(
    icu = 13,
    general = 23,
    hosp = 36,
    deaths_hosp = 5,
    deaths_comm = 3,
    deaths = 8,
    admitted = 53,
    new = 63,
    new_admitted = 116,
    npos_15_64 = 43,
    ntot_15_64 = 83,
    pillar2_pos = 35,
    pillar2_tot = 600,
    pillar2_cases = 35,
    pillar2_over25_pos = 25,
    pillar2_over25_tot = 500,
    pillar2_over25_cases = 25,
    react_pos = 3,
    react_tot = 500)
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
  nms_sero <- c("npos_15_64", "ntot_15_64")
  nms_pillar2 <- c("pillar2_pos", "pillar2_tot")
  nms_pillar2_over25 <- c("pillar2_over25_pos", "pillar2_over25_tot")
  nms_react <- c("react_pos", "react_tot")
  parts <- c(as.list(setdiff(names(observed),
                             c(nms_sero, nms_pillar2,
                               nms_pillar2_over25, nms_react))),
             list(nms_sero), list(nms_pillar2), list(nms_pillar2_over25),
             list(nms_react))

  ll_parts <- lapply(parts, function(x)
    carehomes_compare(state, prev_state, observed_keep(x), pars))

  ## Extremely light testing, though this has already flushed out some
  ## issues
  expect_true(all(lengths(ll_parts) == 6))
  expect_equal(
    carehomes_compare(state, prev_state, observed, pars),
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
  expect_true(all(res[6:17] < population[6:17]))
})


test_that("carehomes_index returns S compartments", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 5, seed = 1L)
  index <- carehomes_index(mod$info())
  expect_equal(
    unname(index$state),
    c(unname(index$run),
      mod$info()$index$hosp_tot,
      mod$info()$index$D_tot,
      mod$info()$index$cum_infections,
      mod$info()$index$S,
      mod$info()$index$cum_admit_by_age))
})


test_that("carehomes_particle_filter_data requires consistent deaths", {
  data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                        1, 0.25)
  ## Add additional columns
  data$deaths_hosp <- data$deaths
  data$deaths_comm <- NA
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
  expect_error(
    carehomes_particle_filter(data),
    "Deaths are not consistently split into total vs community/hospital")
})


test_that("carehomes_particle_filter_data does not allow more than one pillar 2
          data stream", {
  data <- sircovid_data(read_csv(sircovid_file("extdata/example.csv")),
                        1, 0.25)
  ## Add additional columns
  data$deaths_hosp <- NA
  data$deaths_comm <- NA
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

  ## Example that should not cause error
  data$pillar2_pos <- 3
  data$pillar2_tot <- 6
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
})


test_that("the carehomes sircovid model has 19 groups", {
  expect_equal(carehomes_n_groups(), 19)
})


test_that("model_pcr_and_serology_user switch works", {
  set.seed(1)
  p <- carehomes_parameters(0, "england",
                            rel_susceptibility = c(1, 0),
                            vaccine_progression_rate = c(0, 0),
                            waning_rate = 1 / 20,
                            model_pcr_and_serology_user = 0)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  mod$set_state(state)
  mod$set_index(integer(0))

  y <- mod$transform_variables(drop(
    dust::dust_iterate(mod, seq(0, 400, by = 4))))

  ## y$R_neg and y$PCR_neg are increasing over time as noone gets out
  for (i in seq_len(p$n_groups)) {
    for(j in seq_len(p$n_strains))
    {
    expect_true(all(diff(y$R_neg[i, 1, j, ]) >= 0))
    expect_true(all(diff(y$R_neg[i, 2, j, ]) >= 0))
    expect_true(all(diff(y$PCR_neg[i, 1, j, ]) >= 0))
    expect_true(all(diff(y$PCR_neg[i, 2, j, ]) >= 0))
    }
  }

  ## TO DO: ideas for other tests?

})

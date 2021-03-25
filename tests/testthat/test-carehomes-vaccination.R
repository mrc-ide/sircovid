context("carehomes (vaccination)")


test_that("No infections with perfect vaccine wrt rel_susceptibility", {
  ## i.e. if everyone is vaccinated with a vaccine preventing
  ## 100% of acquisition

  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  mod$set_state(state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of s
  expect_equal(length(s), prod(info$dim$S) * 101)
  s <- array(s, c(info$dim$S, 101))

  ## Noone moves into unvaccinated
  ## except in the group where infections because of waning immunity
  expect_true(all(s[-4, 1, ] == 0))

  ## Noone changes compartment within the vaccinated individuals
  expect_true(all(s[, 2, ] == s[, 2, 1]))
})


test_that("No symptomatic infections with perfect vaccine wrt rel_p_sympt", {
  ## i.e. if everyone is vaccinated with a vaccine preventing
  ## 100% of symptoms

  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 0),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  ## Noone moves into I_P, I_C_1 or I_C_2 ever
  ## other than in the 4th age group where some infections are seeded
  ## in the unvaccinated group and because of waning immunity they may
  ## eventually end up in I_P and I_C_2 upon reinfection
  expect_true(all(y$I_P[-4, , , , ] == 0))
  expect_true(all(y$I_C_1[-4, , , , ] == 0))
  expect_true(all(y$I_C_2[-4, , , , ] == 0))

})


test_that("Noone hospitalised with perfect vaccine wrt rel_p_hosp_if_sympt", {
  ## i.e. if everyone is vaccinated with a vaccine preventing
  ## 100% of hospitalisations

  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 0),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  ## Noone moves into hospitalised compartments ever
  ## other than in the 4th age group where some infections are seeded
  ## in the unvaccinated group and because of waning immunity they may
  ## eventually end up in hospital upon reinfection
  expect_true(all(y$D_hosp[-4, ] == 0))
  expect_true(all(y$H_R_unconf[-4, , , , ] == 0))
  expect_true(all(y$H_R_conf[-4, , , , ] == 0))
  expect_true(all(y$H_D_unconf[-4, , , , ] == 0))
  expect_true(all(y$H_D_conf[-4, , , , ] == 0))
})


test_that("No infections with perfect vaccine wrt rel_infectivity", {
  ## i.e. if everyone is vaccinated with a vaccine preventing
  ## 100% of onwards transmission

  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", rel_infectivity = c(1, 0),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  # move susceptibles into vaccinated class
  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  # move initial infections into vaccinated class
  index_I_A <- array(info$index$I_A, info$dim$I_A)
  state[index_I_A[, 1, 1, 2]] <- state[index_I_A[, 1, 1, 1]]
  state[index_I_A[, 1, 1, 1]] <- 0

  mod$set_state(state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of s
  expect_equal(length(s), prod(info$dim$S) * 101)
  s <- array(s, c(info$dim$S, 101))

  ## Noone moves into unvaccinated
  ## except in the group where infections because of waning immunity
  expect_true(all(s[-4, 1, ] == 0))

  ## Noone changes compartment within the vaccinated individuals
  expect_true(all(s[, 2, ] == s[, 2, 1]))
})


test_that("Vaccination of susceptibles works", {
  ## Tests that:
  ## Every susceptible moves to vaccinated and stays there if
  ## everyone quickly gets vaccinated with a vaccine preventing 100% of
  ## acquisition and no waning immunity
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region, beta_value = c(0, 0, 1),
                            beta_date = c(0, 4, 5),
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))
  i <- 4:carehomes_n_groups()

  ## Predefined schedule means we may end up not vaccinating exactly everyone
  ## hence the approximate comparison
  expect_true(all(abs(y$S[i, 1, 1] - y$S[i, 2, 2]) / y$S[i, 1, 1] < 0.01))
  expect_equal(y$S[i, , 101], y$S[i, , 2])
})


test_that("Vaccination of exposed individuals works", {
  ## Tests that:
  ## Every exposed moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if disease progression is stopped after E
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  # stop disease progression after E
  p$gamma_E <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_E <- array(info$index$E, info$dim$E)
  index_S <- array(info$index$S, info$dim$S)
  state[index_E[, , 1, ]] <- round(state[index_S] / 2)
  state[index_E[, , 2, ]] <- round(state[index_S] / 2)
  state[index_S] <- 0

  mod$set_state(state)
  mod$set_index(info$index$E)
  e <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of e
  expect_equal(length(e), prod(info$dim$E) * 101)
  e <- array(e, c(info$dim$E, 101))

  ## every E moves from unvaccinated to vaccinated between
  ## time steps 1 and 2
  E_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  i <- 4:carehomes_n_groups()

  expect_approx_equal(e[i, , E_compartment_idx, unvacc_idx, 1],
                      e[i, , E_compartment_idx, vacc_idx, 2])
  ## then they don't move anymore
  expect_equal(e[i, , E_compartment_idx, vacc_idx, 2],
               e[i, , E_compartment_idx, vacc_idx, 101])
  ## same for second E compartment
  E_compartment_idx <- 2
  expect_approx_equal(e[i, , E_compartment_idx, unvacc_idx, 1],
                      e[i, , E_compartment_idx, vacc_idx, 2])
  ## then they don't move anymore
  expect_equal(e[i, , E_compartment_idx, vacc_idx, 2],
               e[i, , E_compartment_idx, vacc_idx, 101])
})


test_that("Vaccination of asymptomatic infectious individuals works", {
  ## Tests that:
  ## Every I_A moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if disease progression is stopped after I_A
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  # stop disease progression after I_A
  p$gamma_A <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_A[, 1, 1, ]] <- state[index_S]
  state[index_S] <- 0

  mod$set_state(state)
  mod$set_index(info$index$I_A)
  i_A <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of i_A
  expect_equal(length(i_A), prod(info$dim$I_A) * 101)
  i_A <- array(i_A, c(info$dim$I_A, 101))

  ## every I_A moves from unvaccinated to vaccinated between
  ## time steps 1 and 2
  I_A_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  i <- 4:carehomes_n_groups()
  expect_approx_equal(i_A[i, , I_A_compartment_idx, unvacc_idx, 1],
                      i_A[i, , I_A_compartment_idx, vacc_idx, 2])

  ## then they don't move anymore
  expect_equal(
    i_A[i, , I_A_compartment_idx, vacc_idx, 2],
    i_A[i, , I_A_compartment_idx, vacc_idx, 101])
})


test_that("Vaccination of presymptomatic infectious individuals works", {
  ## Tests that:
  ## Every I_P moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if disease progression is stopped after I_P
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  # stop disease progression after I_P
  p$gamma_P <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_I_P <- array(info$index$I_P, info$dim$I_P)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_P[, 1, 1, ]] <- state[index_S]
  state[index_S] <- 0

  mod$set_state(state)
  mod$set_index(info$index$I_P)
  i_P <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of i_P
  expect_equal(length(i_P), prod(info$dim$I_P) * 101)
  i_P <- array(i_P, c(info$dim$I_P, 101))

  ## every I_P moves from unvaccinated to vaccinated between
  ## time steps 1 and 2
  I_P_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  i <- 4:carehomes_n_groups()

  expect_approx_equal(i_P[i, , I_P_compartment_idx, unvacc_idx, 1],
                      i_P[i, , I_P_compartment_idx, vacc_idx, 2])

  ## then they don't move anymore
  expect_equal(
    i_P[i, , I_P_compartment_idx, vacc_idx, 2],
    i_P[i, , I_P_compartment_idx, vacc_idx, 101])
})


test_that("Vaccination of recovered individuals works", {
  ## Test that:
  ## Every R moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if no natural waning of immunity
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_R <- array(info$index$R, info$dim$R)
  index_T_sero_neg <- array(info$index$T_sero_neg, info$dim$T_sero_neg)
  index_T_PCR_neg <- array(info$index$T_PCR_neg, info$dim$T_PCR_neg)
  index_S <- array(info$index$S, info$dim$S)
  state[index_R] <- state[index_S]
  state[index_T_sero_neg] <- state[index_S]
  state[index_T_PCR_neg] <- state[index_S]
  state[index_S] <- 0
  state[index_I_A] <- 0 # remove seeded infections

  mod$set_state(state)
  mod$set_index(info$index$R)
  r <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of r
  expect_equal(length(r), prod(info$dim$R) * 101)
  r <- array(r, c(info$dim$R, 101))

  ## every R moves from unvaccinated to vaccinated between
  ## time steps 1 and 2
  unvacc_idx <- 1
  vacc_idx <- 2
  i <- 4:carehomes_n_groups()
  expect_approx_equal(r[i, , unvacc_idx, 1], r[i, , vacc_idx, 2])
  ## then they don't move anymore
  expect_equal(r[i, , vacc_idx, 2], r[i, , vacc_idx, 101])
})


test_that("Returning to unvaccinated stage works for exposed individuals", {
  ## Tests that:
  ## Every exposed moves back from vaccinated to unvaccinated and stays
  ## there if vaccine has fast waning immunity,
  ## beta is zero, and if disease progression
  ## is stopped after E
  p <- carehomes_parameters(0, "england",
                            beta_value = 0,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, Inf))

  # stop disease progression after E
  p$gamma_E <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_E <- array(info$index$E, info$dim$E)
  index_S <- array(info$index$S, info$dim$S)
  state[index_E[, 1, 1, 2]] <- state[index_S[, 1]]
  state[index_E[, 1, 1, 1]] <- 0
  state[index_E[, 1, 2, 2]] <- state[index_S[, 1]]
  state[index_E[, 1, 2, 1]] <- 0
  state[index_S] <- 0

  mod$set_state(state)
  mod$set_index(info$index$E)
  e <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of e
  expect_equal(length(e), prod(info$dim$E) * 101)
  e <- array(e, c(info$dim$E, 101))

  ## every E moves from vaccinated to unvaccinated between
  ## time steps 1 and 2
  E_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  expect_equal(e[, 1, E_compartment_idx, vacc_idx, 1],
               e[, 1, E_compartment_idx, unvacc_idx, 2])
  ## then they don't move anymore
  expect_equal(e[, 1, E_compartment_idx, unvacc_idx, 2],
               e[, 1, E_compartment_idx, unvacc_idx, 101])
  ## same for second E compartment
  E_compartment_idx <- 2
  expect_equal(e[, 1, E_compartment_idx, vacc_idx, 1],
               e[, 1, E_compartment_idx, unvacc_idx, 2])
  ## then they don't move anymore
  expect_equal(e[, 1, E_compartment_idx, unvacc_idx, 2],
               e[, 1, E_compartment_idx, unvacc_idx, 101])
})


test_that("Returning to unvaccinated stage works for I_A individuals", {
  ## Tests that:
  ## Every I_A moves back from vaccinated to unvaccinated and stays
  ## there if vaccine has fast waning immunity,
  ## beta is zero, and if disease progression
  ## is stopped after I_A
  p <- carehomes_parameters(0, "england",
                            beta_value = 0,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, Inf))

  # stop disease progression after I_A
  p$gamma_A <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_A[, 1, 1, 2]] <- state[index_S[, 1]]
  state[index_I_A[, 1, 1, 1]] <- 0
  state[index_S] <- 0

  mod$set_state(state)
  mod$set_index(info$index$I_A)
  i_A <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of i_A
  expect_equal(length(i_A), prod(info$dim$I_A) * 101)
  i_A <- array(i_A, c(info$dim$I_A, 101))

  ## every I_A moves from vaccinated to unvaccinated between
  ## time steps 1 and 2
  I_A_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  expect_equal(
    i_A[, 1, I_A_compartment_idx, vacc_idx, 1],
    i_A[, 1, I_A_compartment_idx, unvacc_idx, 2])
  ## then they don't move anymore
  expect_equal(
    i_A[, 1, I_A_compartment_idx, unvacc_idx, 2],
    i_A[, 1, I_A_compartment_idx, unvacc_idx, 101])
})


test_that("Returning to unvaccinated stage works for I_P individuals", {
  ## Tests that:
  ## Every I_P moves back from vaccinated to unvaccinated and stays
  ## there if vaccine has fast waning immunity,
  ## beta is zero, and if disease progression
  ## is stopped after I_P
  p <- carehomes_parameters(0, "england",
                            beta_value = 0,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, Inf))

  # stop disease progression after I_P
  p$gamma_P <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_I_P <- array(info$index$I_P, info$dim$I_P)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_P[, 1, 1, 2]] <- state[index_S[, 1]]
  state[index_I_P[, 1, 1, 1]] <- 0
  state[index_S] <- 0

  mod$set_state(state)
  mod$set_index(info$index$I_P)
  i_P <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of i_P
  expect_equal(length(i_P), prod(info$dim$I_P) * 101)
  i_P <- array(i_P, c(info$dim$I_P, 101))

  ## every I_P moves from vaccinated to unvaccinated between
  ## time steps 1 and 2
  I_P_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  expect_equal(
    i_P[, 1, I_P_compartment_idx, vacc_idx, 1],
    i_P[, 1, I_P_compartment_idx, unvacc_idx, 2])
  ## then they don't move anymore
  expect_equal(
    i_P[, 1, I_P_compartment_idx, unvacc_idx, 2],
    i_P[, 1, I_P_compartment_idx, unvacc_idx, 101])
})


test_that("Returning to unvaccinated stage works for recovered individuals", {
  ## Tests that:
  ## Every R moves back from vaccinated to unvaccinated and stays
  ## there if vaccine has fast waning immunity,
  ## beta is zero, and there is no natural immunity
  p <- carehomes_parameters(0, "england",
                            beta_value = 0,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, Inf))

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_R <- array(info$index$R, info$dim$R)
  index_T_sero_neg <- array(info$index$T_sero_neg, info$dim$T_sero_neg)
  index_T_PCR_neg <- array(info$index$T_PCR_neg, info$dim$T_PCR_neg)
  index_S <- array(info$index$S, info$dim$S)
  state[index_R[, 1, 2]] <- state[index_S[, 1]]
  state[index_R[, 1, 1]] <- 0
  state[index_T_sero_neg[, 1, 2]] <- state[index_S[, 1]]
  state[index_T_sero_neg[, 1, 1]] <- 0
  state[index_T_PCR_neg[, 1, 2]] <- state[index_S[, 1]]
  state[index_T_PCR_neg[, 1, 1]] <- 0
  state[index_S] <- 0
  state[index_I_A] <- 0 # remove seeded infections

  mod$set_state(state)
  mod$set_index(info$index$R)
  r <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of r
  expect_equal(length(r), prod(info$dim$R) * 101)
  r <- array(r, c(info$dim$R, 101))

  ## every R moves from vaccinated to unvaccinated between
  ## time steps 1 and 2
  unvacc_idx <- 1
  vacc_idx <- 2
  expect_equal(r[, 1, vacc_idx, 1], r[, 1, unvacc_idx, 2])
  ## then they don't move anymore
  expect_equal(r[, 1, unvacc_idx, 2], r[, 1, unvacc_idx, 101])
})


test_that("Vaccine progression through 3 classes works for susceptibles", {
  ## Tests that:
  ## Every susceptible moves to waning immunity stage and stays there if
  ## everyone quickly gets vaccinated and loses immunity
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            beta_value = 0,
                            rel_susceptibility = c(1, 0, 0),
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(0, Inf, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 3L)

  ## TODO: Anne to look at tidying this parameter up:
  p$model_pcr_and_serology_user <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  i <- 4:carehomes_n_groups()
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))
  expect_approx_equal(y$S[i, 1, 1], y$S[i, 3, 3])
  expect_equal(y$S[i, , 101], y$S[i, , 3])
})


test_that("Vaccine progression through 12 classes works for susceptibles", {
  ## Tests that:
  ## Every susceptible moves to last of 12 waning immunity stage and stays
  ## there if everyone quickly gets vaccinated and loses immunity
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            beta_value = 0,
                            rel_susceptibility = c(1, rep(0, 10), 0),
                            rel_p_sympt = rep(1, 12),
                            rel_p_hosp_if_sympt = rep(1, 12),
                            vaccine_progression_rate =
                              c(0, rep(Inf, 10), 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 12L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 12))))
  i <- 4:carehomes_n_groups()
  expect_approx_equal(y$S[i, 1, 1], y$S[i, 12, 12])
  expect_equal(y$S[i, , 34], y$S[i, , 12])
})


test_that("Clinical progression within a vaccination class works", {
  for (i in 1:50) {
    ## Tests that:
    ## Every susceptible moves to the R compartment corresponding to their
    ## vaccination class if there is a high beta and no vaccine
    ## progression
    p <- carehomes_parameters(0, "england",
                              beta_value = 1e9,
                              rel_susceptibility = c(1, 1, 1),
                              rel_p_sympt = c(1, 1, 1),
                              rel_p_hosp_if_sympt = c(1, 1, 1),
                              vaccine_progression_rate =
                                c(0, 0, 0))

    # increase progression rates
    p[grep("gamma", names(p))] <- 1e9
    # make p_death zero
    p[grep("_D_step", names(p))] <- rep(list(matrix(0)),
                                        length(grep("_D_step", names(p))))

    mod <- carehomes$new(p, 0, 1, seed = 1L)
    info <- mod$info()

    state <- carehomes_initial(info, 1, p)$state

    index_S <- array(info$index$S, info$dim$S)
    index_R <- array(info$index$R, info$dim$R)
    index <- c(index_S, index_R)

    # split S individuals equally between all 3 vaccination groups
    state[index_S[, 2]] <- round(state[index_S[, 1]] / 3)
    state[index_S[, 3]] <- round(state[index_S[, 1]] / 3)
    state[index_S[, 1]] <- state[index_S[, 1]] - state[index_S[, 2]] -
      state[index_S[, 3]]

    mod$set_state(state)
    mod$set_index(index)
    y <- mod$simulate(seq(0, 400, by = 4))

    ## Reshape to show the full shape of s
    expect_equal(length(y),
                 prod(info$dim$S) * 101 + prod(info$dim$R) * 101)
    s <- array(y[seq_len(prod(info$dim$S)), , ], c(info$dim$S, 101))
    r <- array(y[prod(info$dim$S) + seq_len(prod(info$dim$R)), , ],
               c(info$dim$R, 101))

    ## all have moved from S to R in relevant vaccination class
    ## ignoring age group 4 where infections are seeded
    expect_approx_equal(s[-4, , 1], r[-4, 1, , 101])
  }
})


test_that("Returning to unvaccinated stage works for susceptibles", {
  ## Tests that:
  ## Every susceptible moves back to unvaccinated from vacinated if large
  ## waning of immunity and no vaccination
  p <- carehomes_parameters(0, "england",
                            beta_value = 0,
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, Inf))

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  mod$set_state(state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of s
  expect_equal(length(s), prod(info$dim$S) * 101)
  s <- array(s, c(info$dim$S, 101))

  ## noone left in vaccinated at time step 2
  expect_true(all(s[, 2, 2] == 0))

  ## everybody back in unvaccinated at time step 2
  expect_equal(s[, 1, 2], s[, 2, 1])
})


test_that("there are no vaccinated susceptibles when vaccination rate is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  waning_rate <- rep(1 / 20, 19)
  waning_rate[4] <- 0 # no waning in group with seeded infections
  # otherwise S can go up as these infected individuals loose immunity

  p <- carehomes_parameters(0, "england",
                            waning_rate = waning_rate)
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## No vaccinated susceptibles:
  expect_equal(s[-seq_len(carehomes_n_groups()), , ],
               array(0, c(nrow(s) - carehomes_n_groups(), 101)))
})


test_that("Can calculate Rt with an (empty) vaccination class", {
  ## run model with unvaccinated & vaccinated, but both have same susceptibility
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

  ## run model with unvaccinated class only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  rt_1_single_class <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all_single_class <- carehomes_Rt_trajectories(steps, y, p)

  expect_equal(rt_1, rt_1_single_class)
  expect_equal(rt_all, rt_all_single_class)
})


test_that("Effective Rt reduced by rel_susceptibility if all vaccinated", {
  reduced_susceptibility <- 0.2 # can put anything <1 here

  ## run model with unvaccinated & vaccinated (with susceptibility halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, reduced_susceptibility),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

  ## move all individuals to vaccinated
  y_with_vacc <- y
  y_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    y_with_vacc[seq_len(p$n_groups), , ]
  y_with_vacc[seq_len(p$n_groups), , ] <- 0

  rt_1_vacc <- carehomes_Rt(steps, y_with_vacc[, 1, ], p)
  rt_all_vacc <- carehomes_Rt_trajectories(steps, y_with_vacc, p)

  expect_equal(rt_1$eff_Rt_all * reduced_susceptibility,
               rt_1_vacc$eff_Rt_all)
  expect_equal(rt_1$eff_Rt_general * reduced_susceptibility,
               rt_1_vacc$eff_Rt_general)

  expect_equal(rt_all$eff_Rt_all * reduced_susceptibility,
               rt_all_vacc$eff_Rt_all)
  expect_equal(rt_all$eff_Rt_general * reduced_susceptibility,
               rt_all_vacc$eff_Rt_general)

})


test_that("Effective Rt reduced by rel_infectivity if all vaccinated", {
  reduced_infectivity <- 0.2 # can put anything <1 here

  ## run model with unvaccinated & vaccinated (with infectivity halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_infectivity = c(1, reduced_infectivity),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

  ## move all individuals to vaccinated
  y_with_vacc <- y
  y_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    y_with_vacc[seq_len(p$n_groups), , ]
  y_with_vacc[seq_len(p$n_groups), , ] <- 0

  rt_1_vacc <- carehomes_Rt(steps, y_with_vacc[, 1, ], p)
  rt_all_vacc <- carehomes_Rt_trajectories(steps, y_with_vacc, p)

  expect_equal(rt_1$eff_Rt_all * reduced_infectivity,
               rt_1_vacc$eff_Rt_all)
  expect_equal(rt_1$eff_Rt_general * reduced_infectivity,
               rt_1_vacc$eff_Rt_general)

  expect_equal(rt_all$eff_Rt_all * reduced_infectivity,
               rt_all_vacc$eff_Rt_all)
  expect_equal(rt_all$eff_Rt_general * reduced_infectivity,
               rt_all_vacc$eff_Rt_general)

})


test_that("Effective Rt modified if rel_p_sympt is not 1", {
  reduced_p_C <- 0.2 # can put anything <1 here

  ## run model with unvaccinated & vaccinated (with susceptibility halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, reduced_p_C),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)

  ## These are the same as the default values, but setting them again here in
  ## case defaults change as the below assumes mean duration is shorter for
  ## asymptomatic infections
  p$k_A <- 1
  p$gamma_A <- 1 / 2.88
  p$k_P <- 1
  p$gamma_P <- 1 / 1.68
  p$k_C_1 <- 1
  p$gamma_C_1 <- 1 / 2.14
  p$k_C_2 <- 1
  p$gamma_C_2 <- 1 / 1.86

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

  ## move all individuals to vaccinated
  y_with_vacc <- y
  y_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    y_with_vacc[seq_len(p$n_groups), , ]
  y_with_vacc[seq_len(p$n_groups), , ] <- 0

  rt_1_vacc <- carehomes_Rt(steps, y_with_vacc[, 1, ], p)
  rt_all_vacc <- carehomes_Rt_trajectories(steps, y_with_vacc, p)

  # check that the ratio between the Rt with and witout vaccination
  # is constant
  expect_true(all(abs(diff(rt_1_vacc$Rt_all / rt_1$Rt_all)) < 1e-7))

  ## Given mean duration is shorter for asymptomatic individuals, we expect
  ## Rt to be reduced when rel_p_sympt is not 1
  expect_true(all(rt_1_vacc$eff_Rt_all < rt_1$eff_Rt_all))
  expect_true(all(rt_1_vacc$eff_Rt_general < rt_1$eff_Rt_general))

})


test_that("Effective Rt modified if rel_p_hosp_if_sympt is not 1", {
  rel_p_hosp_if_sympt <- 0.2 # can put anything <1 here

  ## run model with unvaccinated & vaccinated (with susceptibility halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, rel_p_hosp_if_sympt),
                            waning_rate = 1 / 20)

  ## set these to non-zero values so that hospitalisation affects R
  p$hosp_transmission <- 0.1
  p$ICU_transmission <- 0.05
  p$G_D_transmission <- 0.05

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

  ## move all individuals to vaccinated
  y_with_vacc <- y
  y_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    y_with_vacc[seq_len(p$n_groups), , ]
  y_with_vacc[seq_len(p$n_groups), , ] <- 0

  rt_1_vacc <- carehomes_Rt(steps, y_with_vacc[, 1, ], p)
  rt_all_vacc <- carehomes_Rt_trajectories(steps, y_with_vacc, p)

  # check that the ratio between the Rt with and witout vaccination
  # is constant
  expect_true(all(abs(diff(rt_1_vacc$Rt_all / rt_1$Rt_all)) < 1e-7))

  ## Given mean duration is shorter for asymptomatic individuals, we expect
  ## Rt to be reduced when rel_p_sympt is not 1
  expect_true(all(rt_1_vacc$eff_Rt_all < rt_1$eff_Rt_all))
  expect_true(all(rt_1_vacc$eff_Rt_general < rt_1$eff_Rt_general))

})


test_that("Can calculate IFR_t with an (empty) vaccination class", {
  ## run model with unvaccinated & vaccinated, but both have same susceptibility
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_S <- mod$info()$index$S
  index_I_weighted <- mod$info()$index$I_weighted
  index <- c(index_S, index_I_weighted)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  S <- y[seq_len(length(index_S)), , ]
  I_weighted <- y[-seq_len(length(index_S)), , ]

  ifr_t_1 <- carehomes_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p)
  ifr_t_all <- carehomes_ifr_t_trajectories(steps, S, I_weighted, p)

  ## run model with unvaccinated class only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_S <- mod$info()$index$S
  index_I_weighted <- mod$info()$index$I_weighted
  index <- c(index_S, index_I_weighted)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  S <- y[seq_len(length(index_S)), , ]
  I_weighted <- y[-seq_len(length(index_S)), , ]

  ifr_t_1_single_class <- carehomes_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p)
  ifr_t_all_single_class <- carehomes_ifr_t_trajectories(steps, S,
                                                         I_weighted, p)

  expect_equal(ifr_t_1, ifr_t_1_single_class)
  expect_equal(ifr_t_all, ifr_t_all_single_class)
})


test_that("IFR_t modified if rel_p_sympt is not 1", {
  reduced_p_C <- 0.2 # can put anything <1 here

  ## run model with unvaccinated & vaccinated (with susceptibility halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, reduced_p_C),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_S <- mod$info()$index$S
  index_I_weighted <- mod$info()$index$I_weighted
  index <- c(index_S, index_I_weighted)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  S <- y[seq_len(length(index_S)), , ]
  I_weighted <- y[-seq_len(length(index_S)), , ]

  ifr_t_1 <- carehomes_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p)
  ifr_t_all <- carehomes_ifr_t_trajectories(steps, S, I_weighted, p)

  ## move all individuals to vaccinated
  S_with_vacc <- S
  S_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    S_with_vacc[seq_len(p$n_groups), , ]
  S_with_vacc[seq_len(p$n_groups), , ] <- 0
  I_weighted_with_vacc <- I_weighted
  I_weighted_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    I_weighted_with_vacc[seq_len(p$n_groups), , ]
  I_weighted_with_vacc[seq_len(p$n_groups), , ] <- 0

  ifr_t_1_vacc <- carehomes_ifr_t(steps, S_with_vacc[, 1, ],
                                  I_weighted_with_vacc[, 1, ], p)
  ifr_t_all_vacc <- carehomes_ifr_t_trajectories(steps, S_with_vacc,
                                                 I_weighted_with_vacc, p)

  ## Given asymptomatic individuals do not die, we expect
  ## IFR_t and IHR_t to be reduced when rel_p_sympt is not 1.
  expect_true(all(ifr_t_1_vacc$IFR_t_all < ifr_t_1$IFR_t_all))
  expect_true(all(ifr_t_1_vacc$IFR_t_general < ifr_t_1$IFR_t_general))
  expect_true(all(ifr_t_1_vacc$IHR_t_all < ifr_t_1$IHR_t_all))
  expect_true(all(ifr_t_1_vacc$IHR_t_general < ifr_t_1$IHR_t_general))

  ## The "no vaccination" version for the model with vaccination should
  ## be the same as the normal version for the model without
  expect_equal(ifr_t_1_vacc$IFR_t_all_no_vacc, ifr_t_1$IFR_t_all)
  expect_equal(ifr_t_1_vacc$IFR_t_general_no_vacc, ifr_t_1$IFR_t_general)
  expect_equal(ifr_t_1_vacc$IHR_t_all_no_vacc, ifr_t_1$IHR_t_all)
  expect_equal(ifr_t_1_vacc$IHR_t_general_no_vacc, ifr_t_1$IHR_t_general)

})


test_that("IFR_t modified if rel_p_hosp_if_sympt is not 1", {
  rel_p_hosp_if_sympt <- 0.2 # can put anything <1 here

  ## run model with unvaccinated & vaccinated
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, rel_p_hosp_if_sympt),
                            waning_rate = 1 / 20)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index_S <- mod$info()$index$S
  index_I_weighted <- mod$info()$index$I_weighted
  index <- c(index_S, index_I_weighted)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(index)
  y <- mod$simulate(steps)
  S <- y[seq_len(length(index_S)), , ]
  I_weighted <- y[-seq_len(length(index_S)), , ]

  ifr_t_1 <- carehomes_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p)
  ifr_t_all <- carehomes_ifr_t_trajectories(steps, S, I_weighted, p)

  ## move all individuals to vaccinated
  S_with_vacc <- S
  S_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    S_with_vacc[seq_len(p$n_groups), , ]
  S_with_vacc[seq_len(p$n_groups), , ] <- 0
  I_weighted_with_vacc <- I_weighted
  I_weighted_with_vacc[seq(p$n_groups + 1, 2 * p$n_groups), , ] <-
    I_weighted_with_vacc[seq_len(p$n_groups), , ]
  I_weighted_with_vacc[seq_len(p$n_groups), , ] <- 0

  ifr_t_1_vacc <- carehomes_ifr_t(steps, S_with_vacc[, 1, ],
                                  I_weighted_with_vacc[, 1, ], p)
  ifr_t_all_vacc <- carehomes_ifr_t_trajectories(steps, S_with_vacc,
                                                 I_weighted_with_vacc, p)

  ## Given asymptomatic individuals do not die, we expect
  ## IFR_t and IHR_t to be reduced when rel_p_sympt is not 1.
  expect_true(all(ifr_t_1_vacc$IFR_t_all < ifr_t_1$IFR_t_all))
  expect_true(all(ifr_t_1_vacc$IFR_t_general < ifr_t_1$IFR_t_general))
  expect_true(all(ifr_t_1_vacc$IHR_t_all < ifr_t_1$IHR_t_all))
  expect_true(all(ifr_t_1_vacc$IHR_t_general < ifr_t_1$IHR_t_general))

  ## The "no vaccination" version for the model with vaccination should
  ## be the same as the normal version for the model without
  expect_equal(ifr_t_1_vacc$IFR_t_all_no_vacc, ifr_t_1$IFR_t_all)
  expect_equal(ifr_t_1_vacc$IFR_t_general_no_vacc, ifr_t_1$IFR_t_general)
  expect_equal(ifr_t_1_vacc$IHR_t_all_no_vacc, ifr_t_1$IHR_t_all)
  expect_equal(ifr_t_1_vacc$IHR_t_general_no_vacc, ifr_t_1$IHR_t_general)

})


test_that("N_tot, N_tot2 and N_tot3 stay constant with vaccination", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  set.seed(1)
  vaccine_schedule <- test_vaccine_schedule(500000, "london")
  p <- carehomes_parameters(0, "london", waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5, 0.1),
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(0, 0, 0.01),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$N_tot3 - mod$transform_variables(y0)$N_tot3 == 0))
  expect_true(all(y$N_tot2 - mod$transform_variables(y0)$N_tot2 == 0))
  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot2 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot3 == 0))
})


test_that("N_tot, N_tot2, N_tot3 are constant with vaccination and k > 1", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  set.seed(1)
  vaccine_schedule <- test_vaccine_schedule(500000, "london")
  p <- carehomes_parameters(0, "london", waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5, 0.1),
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(0, 0, 0.01),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  p[grep("k_", names(p))] <- 2

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$N_tot3 - mod$transform_variables(y0)$N_tot3 == 0))
  expect_true(all(y$N_tot2 - mod$transform_variables(y0)$N_tot2 == 0))
  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot2 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot3 == 0))
})


test_that(
  "N_tot, N_tot2 and N_tot3 stay constant with high rates of vaccination", {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity
    set.seed(1)
    ## TODO: set up a more specific set of tests to test the combined moves
    ## whereby in a single times step an individual progresses to next clinical
    ## stage and progresses to the next vaccination stage
    vaccine_schedule <- test_vaccine_schedule(1000000, "london")
    p <- carehomes_parameters(0, "london", waning_rate = 1 / 20,
                              rel_susceptibility = c(1, 0.5, 0.1),
                              rel_p_sympt = c(1, 1, 1),
                              rel_p_hosp_if_sympt = c(1, 1, 1),
                              vaccine_progression_rate = c(0, 0, 50),
                              vaccine_schedule = vaccine_schedule,
                              vaccine_index_dose2 = 2L)

    mod <- carehomes$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    y0 <- carehomes_initial(info, 1, p)$state
    mod$set_state(carehomes_initial(info, 1, p)$state)
    y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

    expect_true(all(y$N_tot3 - mod$transform_variables(y0)$N_tot3 == 0))
    expect_true(all(y$N_tot2 - mod$transform_variables(y0)$N_tot2 == 0))
    expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
    expect_true(all(colSums(y$N_tot) - y$N_tot2 == 0))
    expect_true(all(colSums(y$N_tot) - y$N_tot3 == 0))
  })

test_that("Outputed vaccination numbers make sense", {
  vaccine_schedule <- test_vaccine_schedule(1000, "london")
  p <- carehomes_parameters(0, "london", waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5, 0.1),
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(0, 0, 0.01),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  ## check outputed objects have correct dimension
  expect_equal(dim(y$cum_n_S_vaccinated), c(19, 3, 101))
  expect_equal(dim(y$cum_n_E_vaccinated), c(19, 3, 101))
  expect_equal(dim(y$cum_n_I_A_vaccinated), c(19, 3, 101))
  expect_equal(dim(y$cum_n_I_P_vaccinated), c(19, 3, 101))
  expect_equal(dim(y$cum_n_R_vaccinated), c(19, 3, 101))

  ## check cumulative stuff is increasing objects have correct dimension

  expect_true(all(apply(y$cum_n_S_vaccinated, c(1, 2), diff) >= 0))
  expect_true(all(apply(y$cum_n_E_vaccinated, c(1, 2), diff) >= 0))
  expect_true(all(apply(y$cum_n_I_A_vaccinated, c(1, 2), diff) >= 0))
  expect_true(all(apply(y$cum_n_I_P_vaccinated, c(1, 2), diff) >= 0))
  expect_true(all(apply(y$cum_n_R_vaccinated, c(1, 2), diff) >= 0))

})

test_that("Outputed S vaccination numbers are what we expect", {
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = c(rep(0, 3), rep(1, 16)))
  p <- carehomes_parameters(0, region, waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  i <- 4:carehomes_n_groups()

  ## there are candidates in S for vaccination
  expect_true(all(y$S[, 1, 1] > 0))
  ## every initial susceptible should be vaccinated within first day
  expect_approx_equal(y$cum_n_S_vaccinated[i, 1, 2], y$S[i, 1, 1])
  ## same for the 10 initially seeded cases
  expect_true(
    all(abs(y$cum_n_I_A_vaccinated[i, 1, 2] - y$I_A[i, 1, 1, 1, 1]) <= 1))

  ## Noone in the first 3 groups vaccinated:
  expect_true(all(y$cum_n_S_vaccinated[1:3, 1, ] == 0))
})


test_that("Outputed E vaccination numbers are what we expect", {
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region, waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  ## stop progression after E to avoid diagonal moves from E to I_A
  p$gamma_E <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  ## empty S ans fill in E initially
  index_E <- array(info$index$E, info$dim$E)
  index_S <- array(info$index$S, info$dim$S)
  state[index_E[, , 1, ]] <- round(state[index_S] / 2)
  state[index_E[, , 2, ]] <- round(state[index_S] / 2)
  state[index_S] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  i <- 4:carehomes_n_groups()
  ## there are candidates in E for vaccination
  expect_true(all(y$E[i, , , 1, 1] > 0))

  ## every initial exposed should be vaccinated within first day
  expect_approx_equal(y$cum_n_E_vaccinated[i, , 2],
                      apply(y$E[i, , , , 1], c(1, 3), sum),
                      rel_tol = 0.15)

  ## same for the 10 initially seeded cases
  expect_true(abs(y$cum_n_I_A_vaccinated[4, 1, 2] -
                  y$I_A[4, , , 1, 1]) <= 3)

})


test_that("Outputed I_A vaccination numbers are what we expect", {
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region, waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  ## stop progression after I_A to avoid diagonal moves from I_A to R
  p$gamma_A <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  ## move people from S to I_A initially
  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_A[, 1, 1, ]] <- state[index_S]
  state[index_S] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  i <- 4:carehomes_n_groups()
  ## there are candidates in I_A for vaccination
  expect_true(all(y$I_A[i, , 1, 1, 1] > 0))
  ## every initial I_A should be vaccinated within first day
  expect_approx_equal(y$cum_n_I_A_vaccinated[i, , 2], y$I_A[i, , , , 1])

})


test_that("Outputed I_P vaccination numbers are what we expect", {
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region, waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  ## stop progression after I_P to avoid diagonal moves out of I_P
  p$gamma_P <- 0

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  ## move people from S to I_P initially
  index_I_P <- array(info$index$I_P, info$dim$I_P)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_P[, 1, , ]] <- state[index_S]
  state[index_S] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  i <- 4:carehomes_n_groups()
  ## there are candidates in I_P for vaccination
  expect_true(all(y$I_P[i, , 1, 1, 1] > 0))
  ## every initial I_P should be vaccinated within first day
  expect_approx_equal(y$cum_n_I_P_vaccinated[i, , 2], y$I_P[i, , , , 1])
  ## same for the 10 initially seeded cases
  expect_equal(y$cum_n_I_A_vaccinated[i, , 2], y$I_A[i, , , , 1])
})


test_that("Outputed R vaccination numbers are what we expect", {
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)
  p <- carehomes_parameters(0, region,
                            waning_rate = 0, # to avoid diagonal moves out of R
                            rel_susceptibility = c(1, 0.5),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  ## empty S ans fill in R initially
  index_R <- array(info$index$R, info$dim$R)
  index_T_sero_neg <- array(info$index$T_sero_neg, info$dim$T_sero_neg)
  index_T_PCR_neg <- array(info$index$T_PCR_neg, info$dim$T_PCR_neg)
  index_S <- array(info$index$S, info$dim$S)
  state[index_R[, , ]] <- state[index_S]
  state[index_T_sero_neg[, , ]] <- state[index_S]
  state[index_T_PCR_neg[, , ]] <- state[index_S]
  state[index_S] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))

  i <- 4:carehomes_n_groups()
  ## there are candidates in R for vaccination
  expect_true(all(y$R[i, , 1, 1] > 0))
  ## every initial recovered should be vaccinated within first day
  expect_approx_equal(y$cum_n_R_vaccinated[i, , 2], y$R[i, , , 1])
  ## same for the 10 initially seeded cases
  expect_equal(y$cum_n_I_A_vaccinated[i, , 2], y$I_A[i, , , , 1])
})


test_that("check_rel_param rejects out of bounds errors", {
  expect_error(
    check_rel_param(NULL, "rel_param"),
    "At least one value required for rel_param")
  expect_error(
    check_rel_param(-1, "rel_param"),
    "All values of rel_param must lie in [0, 1]",
    fixed = TRUE)
  expect_error(
    check_rel_param(1.1, "rel_param"),
    "All values of rel_param must lie in [0, 1]",
    fixed = TRUE)
  expect_error(
    check_rel_param(c(1, 1.8), "rel_param"),
    "All values of rel_param must lie in [0, 1]",
    fixed = TRUE)
  expect_error(
    check_rel_param(t(c(0.9, 0.8)), "rel_param"),
    "First value of rel_param must be 1")
})


test_that("check_rel_param allows sensible inputs", {
  expect_silent(
    check_rel_param(t(c(1, 0.5, 0.7)), "rel_param"))
  expect_silent(
    check_rel_param(t(c(1, 0.7, 0.5)), "rel_param"))
  expect_silent(
    check_rel_param(t(1), "rel_param"))
  expect_silent(
    check_rel_param(t(c(1, 1)), "rel_param"))
  expect_silent(
    check_rel_param(t(c(1, 0)), "rel_param"))
  expect_silent(
    check_rel_param(t(c(1, 0, 1)), "rel_param"))
})


test_that("build_rel_param rejects wrong dimension or out of bound inputs", {
  expect_error(
    build_rel_param(
      rel_param = matrix(c(1, 0.5, 1, 0.7), nrow = 2, byrow = TRUE),
      n_vacc_classes = 2, "rel_param"),
    "rel_param should have as many rows as age groups")
  expect_error(
    build_rel_param(10, n_vacc_classes = 1, "rel_param"),
    "All values of rel_param must lie in [0, 1]", fixed = TRUE)
})


test_that("build_rel_param works as expected", {
  expect_equal(
    build_rel_param(1, n_vacc_classes = 1, "rel_param"),
    matrix(1, nrow = carehomes_n_groups(), ncol = 1))
  mat <- matrix(rep(c(1, 0.1), carehomes_n_groups()), byrow = TRUE,
                nrow = carehomes_n_groups(), ncol = 2)
  expect_equal(
    build_rel_param(c(1, 0.1), n_vacc_classes = 2, "rel_param"),
    mat)
  expect_equal(
    build_rel_param(mat, n_vacc_classes = 2, "rel_param"),
    mat)
  mat_rand <- cbind(rep(1, carehomes_n_groups()), runif(carehomes_n_groups()))
  expect_equal(
    build_rel_param(mat_rand, n_vacc_classes = 2, "rel_param"),
    mat_rand)
})


test_that("build_vaccine_progression_rate rejects insensible inputs", {
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate =
                                     cbind(rep(-1, 19), rep(1, 19), rep(1, 19)),
                                   n_vacc_classes = 3),
    "'vaccine_progression_rate' must have only non-negative values")
  msg1 <- "'vaccine_progression_rate' must be either:"
  msg2 <- "a vector of length 'n_vacc_classes'"
  msg3 <- "or a matrix with 'n_groups' rows and 'n_vacc_classes' columns"
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = c(1, 1, 1),
                                   n_vacc_classes = 2),
    paste(msg1, msg2, msg3))
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = c(1, 1, -1),
                                   n_vacc_classes = 3),
    "'vaccine_progression_rate' must have only non-negative values")
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = matrix(1, 19, 5),
                                   n_vacc_classes = 3),
    "'vaccine_progression_rate' must have 'n_vacc_classes' columns")
  expect_error(
    build_vaccine_progression_rate(vaccine_progression_rate = matrix(1, 9, 3),
                                   n_vacc_classes = 3),
    "'vaccine_progression_rate' must have as many rows as age groups")
  ntot <- rep(1000, 19)
  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      rel_susceptibility = c(1, 1, 1),
      vaccine_progression_rate = c(1, 1, 1)),
    "Column 1 of 'vaccine_progression_rate' must be zero (dose 1)",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      rel_susceptibility = c(1, 1, 1),
      vaccine_progression_rate = matrix(c(1, 1, 1), 19, 3)),
    "Column 1 of 'vaccine_progression_rate' must be zero (dose 1)",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      rel_susceptibility = 1,
      vaccine_progression_rate = 1),
    "Column 1 of 'vaccine_progression_rate' must be zero (dose 1)",
    fixed = TRUE)

  schedule <- vaccine_schedule(date = 1L,
                               doses = array(0, c(19, 2, 5)))
  dt <- 1 / 4
  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      dt,
      rel_susceptibility = c(1, 1, 1),
      vaccine_progression_rate = c(0, 1, 1),
      vaccine_index_dose2 = 2,
      vaccine_schedule = schedule),
    "Column 2 of 'vaccine_progression_rate' must be zero (dose 2)",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      dt,
      rel_susceptibility = c(1, 1, 1, 1),
      vaccine_progression_rate = c(0, 1, 1, 1),
      vaccine_index_dose2 = 3,
      vaccine_schedule = schedule),
    "Column 3 of 'vaccine_progression_rate' must be zero (dose 2)",
    fixed = TRUE)
})


test_that("build_vaccine_progression_rate allows sensible inputs and works", {
  expect_silent(
    build_vaccine_progression_rate(vaccine_progression_rate = 0,
                                   n_vacc_classes = 1,
                                   index_dose = c(1, 1)))
  expect_equal(
    build_vaccine_progression_rate(vaccine_progression_rate = c(0, 1),
                                   n_vacc_classes = 2,
                                   index_dose = c(1, 1)),
    cbind(rep(0, 19), rep(1, 19)))
  expect_silent(
    build_vaccine_progression_rate(vaccine_progression_rate = matrix(0, 19, 3),
                                   n_vacc_classes = 3,
                                   index_dose = c(1, 1)))
  expect_silent(
    build_vaccine_progression_rate(vaccine_progression_rate = NULL,
                                   n_vacc_classes = 3,
                                   index_dose = c(1, 1)))
  expect_equal(
    build_vaccine_progression_rate(vaccine_progression_rate = NULL,
                                   n_vacc_classes = 3,
                                   index_dose = c(1, 1)),
    matrix(0, 19, 3))
})


test_that("build_waning_rate works as expected", {
  expect_error(
    build_waning_rate(NULL),
    "At least one value required for 'waning_rate'")
  expect_error(
    build_waning_rate(-1),
    "'waning_rate' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    build_waning_rate(rep(0.5, 2)),
    "'waning_rate' should have as many elements as age groups")
  expect_equal(
    build_waning_rate(0),
    rep(0, 19))
  expect_equal(
    build_waning_rate(0.5),
    rep(0.5, 19))
  expect_equal(
    build_waning_rate(rep(0.5, 19)),
    rep(0.5, 19))
})


## Heading towards real-life use, let's vaccinate people at a rate of
## 5k/day. This does not run an epidemic beforehand though, and we'll
## use a "null" vaccine for now.
test_that("run sensible vaccination schedule", {
  region <- "east_of_england"
  uptake <- c(rep(0, 3), rep(1, 16))
  daily_doses <- rep(50000, 120)
  n <- vaccine_priority_population(region, uptake,
                                   prop_hcw = rep(0, 19),
                                   prop_very_vulnerable = rep(0, 19),
                                   prop_underlying_condition = rep(0, 19))
  vaccine_schedule <- vaccine_schedule_future(0, daily_doses, 200, n)
  expect_equal(sum(vaccine_schedule$doses[, 2, ]), 0)
  expect_equal(sum(vaccine_schedule$doses[1:3, , ]), 0)

  p <- carehomes_parameters(0, "east_of_england",
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)
  ## TODO: Anne to look at tidying this parameter up:
  p$model_pcr_and_serology_user <- 0

  ## Let's go:
  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state
  ## Remove seed, so that we have no infection process here:
  state[state == 10] <- 0

  mod$set_state(state)
  mod$set_index(integer(0))

  keep <- c("cum_n_S_vaccinated",
            "cum_n_E_vaccinated",
            "cum_n_I_A_vaccinated",
            "cum_n_I_P_vaccinated",
            "cum_n_R_vaccinated")
  index <- unlist(lapply(info$index[keep], "[", 1:19), FALSE, FALSE)

  mod$set_index(index)
  y <- mod$simulate(seq(0, 380, by = 4)[-1])
  s <- array(y, c(19, 5, dim(y)[3]))

  ## Never vaccinate any young person:
  expect_true(all(s[1:3, , ] == 0))

  ## Sum over compartments
  cum_n_vaccinated <- t(apply(s, c(1, 3), sum))
  n_vaccinated <- diff(cum_n_vaccinated)

  ## You can visualise the vaccination process here:
  ## > matplot(m, type = "l", lty = 1)

  tot <- rowSums(n_vaccinated)
  expect_true(all(tot >= 49000 & tot < 51000))

  ## Vaccinate all the CHW/CHR first, then down the priority
  ## groups. This is easy to check visually but harder to describe:
  priority <- list(18:19, 17, 16, 15, 14, 13, 12, 11,
                   9:10, 7:8, 1:6)
  i <- lapply(priority, function(p)
    range(c(apply((n_vaccinated > 5000)[, p, drop = FALSE], 2, which))))
  for (j in seq_along(i)) {
    if (j > 2) {
      ## using <= as if many doses available each day you may vaccinate
      ## several priority groups in the same day
      expect_true(max(unlist(i[seq_len(j - 2)])) <= i[[j]][[1]])
    }
    if (j > 1) {
      expect_true(all(i[[j - 1]][[1]] <= i[[j]][[1]]))
    }
  }
})


test_that("can add vaccination to a set of model state", {
  region <- "east_of_england"

  p_orig <- carehomes_parameters(0, region)

  vaccine_schedule <- test_vaccine_schedule(daily_doses = 5000,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)

  p_vacc <- carehomes_parameters(0, region,
                                 rel_susceptibility = c(1, 1),
                                 rel_p_sympt = c(1, 1),
                                 rel_p_hosp_if_sympt = c(1, 1),
                                 vaccine_schedule = vaccine_schedule,
                                 vaccine_index_dose2 = 2L)

  mod_orig <- carehomes$new(p_orig, 0, 10, seed = 1L)
  mod_vacc <- carehomes$new(p_vacc, 0, 10, seed = 1L)
  state_orig <- mod_orig$state()
  state_orig[] <- seq_along(state_orig)

  state_vacc <- vaccine_remap_state(state_orig, mod_orig$info(),
                                        mod_vacc$info())
  expect_equal(sum(state_vacc), sum(state_orig))

  tmp_orig <- mod_orig$transform_variables(state_orig)
  tmp_vacc <- mod_vacc$transform_variables(state_vacc)
  cmp <- function(v_orig, v_vacc, name) {
    nd <- length(dim(v_orig))
    if (identical(dim(v_orig), dim(v_vacc))) {
      identical(v_orig, v_vacc)
    } else if (nd == 2) {
      identical(v_orig, v_vacc[, 1, drop = FALSE]) &&
        all(v_vacc[, -1] == 0)
    } else if (nd == 3) {
      identical(v_orig, v_vacc[, 1, , drop = FALSE]) &&
        all(v_vacc[, -1, ] == 0)
    } else if (nd == 4) {
      identical(v_orig, v_vacc[, , 1, , drop = FALSE]) &&
        all(v_vacc[, , -1, ] == 0)
    }
  }
  expect_true(all(unlist(Map(cmp, tmp_orig, tmp_vacc, names(tmp_orig)))))
})


test_that("vaccine_uptake must be the correct length", {
  expect_error(
    vaccine_priority_proportion(uptake = c(0, 0, 0)),
    "Invalid length 3 for 'uptake', must be 1 or 19")
})


test_that("If vaccine dose 2 index provided, we need a schedule", {
  schedule <- vaccine_schedule(date = 1L,
                               doses = array(0, c(19, 2, 5)))
  dt <- 1 / 4
  ntot <- rep(1000, 19)
  expect_error(
    carehomes_parameters_vaccination(ntot, vaccine_index_dose2 = 2L),
    "'vaccine_index_dose2' set without schedule")
})


test_that("vaccine_index_dose2 must have valid value", {
  schedule <- vaccine_schedule(date = 1L,
                               doses = array(0, c(19, 2, 5)))
  dt <- 1 / 4
  ntot <- rep(1000, 19)
  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      dt,
      vaccine_index_dose2 = 2L,
      vaccine_schedule = schedule),
    "Invalid value for 'vaccine_index_dose2', must be in [1, 1]",
    fixed = TRUE)

  expect_error(
    carehomes_parameters_vaccination(
      ntot,
      dt,
      rel_susceptibility = c(1, 1, 1, 1),
      vaccine_progression_rate = c(0, 1, 1, 1),
      vaccine_index_dose2 = 6L,
      vaccine_schedule = schedule),
    "Invalid value for 'vaccine_index_dose2', must be in [1, 4]",
    fixed = TRUE)
})


test_that("can upgrade model state", {
  ## This is the situation upgrading sircovid 0.7.2 -> 0.8.0 as we
  ## lack cum_n_vaccinated; that can obviously be added.

  region <- "east_of_england"

  p_orig <- carehomes_parameters(0, region)
  mod_orig <- carehomes$new(p_orig, 0, 10, seed = 1L)
  info_orig <- mod_orig$info()

  ## Elimate our index
  info_orig$index <-
    info_orig$index[names(info_orig$index) != "cum_n_vaccinated"]
  info_orig$dim <-
    info_orig$dim[names(info_orig$dim) != "cum_n_vaccinated"]
  ## Reset the index
  k <- 0L
  for (i in seq_along(info_orig$index)) {
    n <- length(info_orig$index[[i]])
    info_orig$index[[i]] <- seq_len(n) + k
    k <- k + n
  }
  info_orig$len <- k

  state_orig <- matrix(0, info_orig$len, 10)
  state_orig[] <- seq_along(state_orig)

  ## Then a new set with vaccination:
  vaccine_schedule <- test_vaccine_schedule(daily_doses = 5000,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)

  p_vacc <- carehomes_parameters(0, region,
                                 rel_susceptibility = c(1, 1),
                                 rel_p_sympt = c(1, 1),
                                 rel_p_hosp_if_sympt = c(1, 1),
                                 vaccine_schedule = vaccine_schedule,
                                 vaccine_index_dose2 = 2L)

  mod_vacc <- carehomes$new(p_vacc, 0, 10, seed = 1L)
  info_vacc <- mod_vacc$info()
  state_vacc <- vaccine_remap_state(state_orig, info_orig, info_vacc)

  expect_equal(sum(state_vacc), sum(state_orig))
  expect_equal(state_vacc[info_vacc$index$cum_n_vaccinated],
               rep(0, 2 * carehomes_n_groups()))
})


test_that("Refuse to upgrade impossible model state", {
  region <- "east_of_england"

  p_orig <- carehomes_parameters(0, region)
  mod_orig <- carehomes$new(p_orig, 0, 10, seed = 1L)
  info_orig <- mod_orig$info()

  ## Elimate our index
  info_orig$index <-
    info_orig$index[names(info_orig$index) != "E"]
  info_orig$dim <-
    info_orig$dim[names(info_orig$dim) != "E"]
  ## Reset the index
  k <- 0L
  for (i in seq_along(info_orig$index)) {
    n <- length(info_orig$index[[i]])
    info_orig$index[[i]] <- seq_len(n) + k
    k <- k + n
  }
  info_orig$len <- k

  state_orig <- matrix(0, info_orig$len, 10)
  state_orig[] <- seq_along(state_orig)

  ## Then a new set with vaccination:
  vaccine_schedule <- test_vaccine_schedule(daily_doses = 5000,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)

  p_vacc <- carehomes_parameters(0, region,
                                 rel_susceptibility = c(1, 1),
                                 rel_p_sympt = c(1, 1),
                                 rel_p_hosp_if_sympt = c(1, 1),
                                 vaccine_schedule = vaccine_schedule,
                                 vaccine_index_dose2 = 2L)

  mod_vacc <- carehomes$new(p_vacc, 0, 10, seed = 1L)
  info_vacc <- mod_vacc$info()
  expect_error(
    vaccine_remap_state(state_orig, info_orig, info_vacc),
    "Can't remap state (can't add variables 'E')",
    fixed = TRUE)
  expect_error(
    vaccine_remap_state(matrix(0, info_vacc$len, 10), info_vacc, info_orig),
    "Can't downgrade state (previously had variables 'E')",
    fixed = TRUE)
})


test_that("Can vaccinate given a schedule", {
  p <- carehomes_parameters(0, "england", rel_susceptibility = c(1, 1, 0),
                            beta_value = 0,
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(0, 0, 0),
                            waning_rate = 1 / 20)
  p$index_dose <- c(1L, 2L)
  end_date <- sircovid_date("2020-06-01")

  start_vacc_date_1 <- sircovid_date("2020-03-01")
  delay_vacc_date_2 <- 28
  ndays_vacc <- 31
  i <- seq(start_vacc_date_1 * 4, length.out = ndays_vacc * 4)
  m <- array(0, c(19, 2, (end_date + 1) * 4))
  step_doses_1 <- 10000
  step_doses_2 <- 2000
  m[, 1, i] <- step_doses_1 # first dose schedule
  m[, 2, i + delay_vacc_date_2 * 4] <- step_doses_2 # second dose schedule
  p$vaccine_dose_step <- m

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  mod$set_state(state)
  steps <- seq(0, end_date * 4, by = 4)
  y <- mod$transform_variables(mod$simulate(steps))

  #### check first dose schedule
  n_vacc_fisrt_dose <- apply(y$cum_n_vaccinated[, 1, 1, ], 1, diff)

  ## check no vaccination before wanted date
  days_vacc_1 <- seq(start_vacc_date_1 + 1, start_vacc_date_1 + ndays_vacc)
  expect_true(
    all(n_vacc_fisrt_dose[seq_len(days_vacc_1[1] - 1), ] == 0))

  ## check that in all groups but 18:19
  ## (which are small and therefore get vaccinated faster)
  ## we get the right number of vaccinations per day in the wanted interval
  x <- n_vacc_fisrt_dose[days_vacc_1, - (18:19)]
  expect_approx_equal(x, matrix(step_doses_1 * 4, nrow(x), ncol(x)))

  ## check that no vaccination after wanted date
  ## we get the right number of vaccinations per day in the wanted interval
  expect_true(
    all(n_vacc_fisrt_dose[seq(last(days_vacc_1) + 1, end_date, 1), ] == 0))

  #### check second dose schedule
  n_vacc_second_dose <- apply(y$cum_n_vaccinated[, 2, 1, ], 1, diff)

  ## check no vaccination before wanted date
  days_vacc_2 <- seq(start_vacc_date_1 + delay_vacc_date_2 + 1,
                     start_vacc_date_1 + delay_vacc_date_2 + ndays_vacc)
  expect_true(
    all(n_vacc_second_dose[seq_len(days_vacc_2[1] - 1), ] == 0))

  ## check that
  ## we get the right number of vaccinations per day in the wanted interval
  x <- n_vacc_second_dose[days_vacc_2, ]
  expect_approx_equal(x, matrix(step_doses_2 * 4, nrow(x), ncol(x)))

  ## check that no vaccination after wanted date
  ## we get the right number of vaccinations per day in the wanted interval
  expect_true(
    all(n_vacc_second_dose[seq(last(days_vacc_2) + 1, end_date, 1), ] == 0))

})


test_that("can create parameters with vaccination data", {
  region <- "london"
  uptake_by_age <- test_example_uptake()
  daily_doses <- rep(20000, 365)
  mean_days_between_doses <- 12 * 7
  n <- vaccine_priority_population(region, uptake_by_age)
  date_start_vaccination <- sircovid_date("2020-02-01")
  schedule <- vaccine_schedule_future(
    date_start_vaccination, daily_doses, mean_days_between_doses, n)

  p <- carehomes_parameters(0, region,
                            rel_susceptibility = c(1, 1, 0),
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            waning_rate = 1 / 20,
                            vaccine_index_dose2 = 2L,
                            vaccine_schedule = schedule)
  expect_equal(dim(p$vaccine_dose_step),
               c(19, 2, (date_start_vaccination + length(daily_doses)) * 4))
})

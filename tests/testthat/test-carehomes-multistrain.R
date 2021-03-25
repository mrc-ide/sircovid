context("carehomes (multistrain)")

test_that("carehomes_parameters_strain works as expected", {
  expect_error(
    carehomes_parameters_strain(NULL, NULL, NULL, 1),
    "At least one value required for 'strain_transmission'")
  expect_error(
    carehomes_parameters_strain(c(1, -1), NULL, NULL, 1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 2), NULL, NULL, 1),
    "'strain_transmission[[1]]' must be 1",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 1), NULL, NULL, 1),
    "'strain_transmission[[1]]' must be 1",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 3), NULL, NULL, 1),
    "Only 1 or 2 strains valid",
    fixed = TRUE)
  expect_equal(
    carehomes_parameters_strain(1, NULL, NULL, 1),
    list(n_strains = 1,
         strain_transmission = 1,
         strain_seed_step = 0))
  expect_equal(
    carehomes_parameters_strain(c(1, 1), NULL, NULL, 1),
    list(n_strains = 2,
         strain_transmission = c(1, 1),
         strain_seed_step = 0))
  expect_equal(
    carehomes_parameters_strain(c(1, 2), NULL, NULL, 1),
    list(n_strains = 2,
         strain_transmission = c(1, 2),
         strain_seed_step = 0))
})


test_that("Prevent impossible seedings", {
  expect_error(
    carehomes_parameters_strain(c(1, 1), NULL, 1, 0.25),
    "As 'strain_seed_date' is NULL, expected 'strain_seed_rate' to be NULL")
  expect_error(
    carehomes_parameters_strain(1, c(10, 20), 1, 0.25),
    "Can't use 'strain_seed_date' if only using one strain")
  expect_error(
    carehomes_parameters_strain(c(1, 1), c(10, 20), c(-1, 0), 1),
    "'strain_seed_rate' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(c(1, 1), c(10, 20, 30), 1, 0.25),
    "'strain_seed_date' and 'strain_seed_rate' must be the same length")
})


test_that("Can seed with one-day window", {
  date <- c("2020-03-01", "2020-03-02")
  rate <- c(100, 0)
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), rate, 1 / 4)
  expect_equal(sum(p$strain_seed_step), 100)
  expect_equal(tail(p$strain_seed_step, 6), c(0, 25, 25, 25, 25, 0))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-02"))
})


test_that("Can seed with ten-day window", {
  date <- c("2020-03-01", "2020-03-11")
  rate <- c(100, 0)
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), rate, 1 / 4)
  expect_equal(sum(p$strain_seed_step), 100 * 10)
  expect_equal(tail(p$strain_seed_step, 6), c(25, 25, 25, 25, 25, 0))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-11"))
})


test_that("Can seed with > 2 dates", {
  date <- c("2020-03-01", "2020-03-10", "2020-03-20", "2020-03-21")
  rate <- c(100, 5, 20, 1)
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), rate, 1 / 4)
  expect_equal(as.numeric(table(p$strain_seed_step)),
               c(sircovid_date("2020-03-01") * 4 - 1, 1,  40, 4, 36))
  expect_equal(tail(p$strain_seed_step, 6), c(1.25, 5, 5, 5, 5, 0.25))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-21"))
})


test_that("Adding empty strains makes no difference", {
  p1 <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 0))
  p3 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 0.5))
  np <- 10
  mod1 <- carehomes$new(p1, 0, np, seed = 1L)
  mod2 <- carehomes$new(p2, 0, np, seed = 1L)
  mod3 <- carehomes$new(p3, 0, np, seed = 1L)
  end <- sircovid_date("2020-03-31") / p1$dt

  mod1$set_index(carehomes_index(mod1$info())$run)
  mod2$set_index(carehomes_index(mod2$info())$run)
  mod3$set_index(carehomes_index(mod3$info())$run)

  initial1 <- carehomes_initial(mod1$info(), 1, p1)
  initial2 <- carehomes_initial(mod1$info(), 1, p2)
  initial3 <- carehomes_initial(mod1$info(), 1, p3)

  res1 <- mod1$run(end)
  res2 <- mod2$run(end)
  res3 <- mod3$run(end)

  expect_equal(res2, res1)
  expect_equal(res3, res1)
})


test_that("Seeding of second strain generates an epidemic", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  ## Did the seeded cases go on to infect other people?
  expect_true(y$cum_infections_per_strain[2, 101] > n_seeded_new_strain_inf)

  ## Check the epidemic of the second strain starts when we expect
  steps <- seq(0, 400, by = 4)
  date <- sircovid_date_as_date(steps / 4)
  s_date <- sircovid_date(date)
  s_date_seeding <- sircovid_date(date_seeding)

  ## No cases before seeding
  expect_true(all(y$E[, 2, , , s_date < s_date_seeding] == 0))
  ## No cases on seeding day other than in 4th age group
  expect_true(all(y$E[-4, 2, , , s_date == s_date_seeding] == 0))
  ## Some cases on seeding day in 4th age group
  expect_true(y$E[4, 2, 1, , s_date == s_date_seeding] > 0)
  ## Some cases on all days after seeding day

  ## It's not guaranteed that *all* will be greater than zero, but most will be
  expect_true(mean(colSums(y$E[, 2, 1, , s_date >= s_date_seeding]) > 0) > 0.9)
})


test_that("Second more virulent strain takes over", {
  np <- 10
  n_seeded_new_strain_inf <- 10
  start_date <- sircovid_date("2020-02-07")
  date_seeding <- start_date # seed both strains on same day
  p <- carehomes_parameters(start_date, "england",
                            strain_transmission = c(1, 10),
                            strain_seed_date = c(date_seeding,
                                                 date_seeding + 1),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  mod <- carehomes$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## cumulative infections with 2nd strain larger than with 1st strain
  ## (average over 10 runs)
  expect_true(mean(y$cum_infections_per_strain[1, , 101]) <
                mean(y$cum_infections_per_strain[2, , 101]))
})


test_that("Second less virulent strain does not take over", {
  np <- 10
  n_seeded_new_strain_inf <- 10
  start_date <- sircovid_date("2020-02-07")
  date_seeding <- start_date # seed both strains on same day
  p <- carehomes_parameters(start_date, "england",
                            strain_transmission = c(1, 0.1),
                            strain_seed_date = c(date_seeding,
                                                 date_seeding + 1),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  mod <- carehomes$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  ## Cumulative infections with 2nd strain smaller than with 1st strain
  ## (average over 10 runs)
  expect_true(mean(y$cum_infections_per_strain[1, , 101]) >
                mean(y$cum_infections_per_strain[2, , 101]))

})


test_that("N_tot, N_tot2 and N_tot3 stay constant with second strain", {
  ## Default for waning_rate is 0, setting to a non-zero value so that
  ## this test passes with waning immunity
  set.seed(1)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            waning_rate = 1 / 20,
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$N_tot3 - mod$transform_variables(y0)$N_tot3 == 0))
  expect_true(all(y$N_tot2 - mod$transform_variables(y0)$N_tot2 == 0))
  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot2 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot3 == 0))
})


test_that("No infection after seeding of second strain with 0 transmission", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## can't compare to fixed number so generate a feasible range
  pois_range <- rpois(1000, n_seeded_new_strain_inf)

  ## Expect the seeded cases did not infect any other people
  expect_true(y$cum_infections_per_strain[2, 101] <= max(pois_range))
  expect_true(y$cum_infections_per_strain[2, 101] >= min(pois_range))
})


test_that("Everyone is infected when second strain transmission is large", {
  n_seeded_new_strain_inf <- 10
  date_seeding <- c("2020-03-07", "2020-03-08")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1e9),
                            strain_seed_date =
                              sircovid_date(date_seeding),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  ## set gamma_E to Inf so that seeded individuals move through each E stage
  ## in one step
  p$gamma_E <- Inf

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  steps <- seq(0, 400, by = 4)
  date <- sircovid_date_as_date(steps / 4)
  s_date <- sircovid_date(date)
  s_date_seeding <- sircovid_date(date_seeding[[1]])
  ## No cases before seeding
  expect_true(all(y$E[, 2, , , s_date < s_date_seeding] == 0))

  ## The +2 is because we need seeded individuals to get out of the first and
  ## second E compartments before they can go on to infect others
  expect_true(all(y$S[, 1, s_date > (s_date_seeding + 2)] == 0))
})


test_that("No infection with either strain with perfect vaccine", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england",
                            strain_transmission = c(1, 1),
                            rel_susceptibility = c(1, 0),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  index_E <- array(info$index$E, info$dim$E)
  state[index_E[4, 2, 1, 1]] <- 10 # seed infections with second strain

  mod$set_state(state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## Noone moves into unvaccinated
  ## except in the group where infections because of waning immunity
  expect_true(all(y$S[-4, 1, ] == 0))

  ## Noone changes compartment within the vaccinated individuals
  expect_true(all(y$S[, 2, ] == y$S[, 2, 1]))

  ## Noone gets infected with either strain
  expect_true(all(y$cum_infections_per_strain == 0))
})


test_that("different strains are equivalent", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1))
  np <- 10
  mod <- carehomes$new(p, 0, np, seed = 1L, n_threads = 10)
  end <- sircovid_date("2020-07-31") / p$dt
  end <- sircovid_date("2020-03-31") / p$dt

  initial <- carehomes_initial(mod$info(), 1, p)
  y <- mod$transform_variables(initial$state)
  y$I_A <- y$I_A[, 2:1, , , drop = FALSE]
  y$T_PCR_pos <- y$T_PCR_pos[, 2:1, , , drop = FALSE]
  y$T_sero_pre <- y$T_sero_pre[, 2:1, , , drop = FALSE]

  initial2_state <- unlist(y)

  mod$set_state(initial$state, initial$step)
  index <- mod$info()$index
  index_run <- c(icu = index[["ICU_tot"]],
                 general = index[["general_tot"]],
                 deaths_comm = index[["D_comm_tot"]],
                 deaths_hosp = index[["D_hosp_tot"]],
                 admitted = index[["cum_admit_conf"]],
                 new = index[["cum_new_conf"]],
                 sero_pos = index[["sero_pos"]],
                 sympt_cases = index[["cum_sympt_cases"]],
                 sympt_cases_over25 = index[["cum_sympt_cases_over25"]],
                 react_pos = index[["react_pos"]],
                 infections = index[["cum_infections"]])

  steps <- seq(initial$step, end, by = 4)
  mod$set_index(index_run)
  res1 <- mod$simulate(steps)

  mod2 <- carehomes$new(p, 0, np, seed = 1L, n_threads = 10)
  mod2$set_state(initial2_state, initial$step)
  mod2$set_index(index_run)
  res2 <- mod2$simulate(steps)

  expect_equal(res1, res2)
})


test_that("Swapping strains gives identical results with different index", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1))

  np <- 1
  mod <- carehomes$new(p, 0, np, seed = 1L)
  end <- sircovid_date("2020-05-1") / p$dt
  initial <- carehomes_initial(mod$info(), 1, p)
  y <- mod$transform_variables(initial$state)
  y$I_A <- y$I_A[, 2:1, , , drop = FALSE]
  y$T_PCR_pos <- y$T_PCR_pos[, 2:1, , , drop = FALSE]
  y$T_sero_pre <- y$T_sero_pre[, 2:1, , , drop = FALSE]

  initial2_state <- unlist(y)
  mod$set_state(initial$state, initial$step)
  index <- mod$info()$index

  steps <- seq(initial$step, end, by = 1)

  res1 <- drop(mod$simulate(steps))

  mod2 <- carehomes$new(p, 0, np, seed = 1L)

  mod2$set_state(initial2_state, initial$step)
  res2 <- drop(mod2$simulate(steps))

  z1 <- mod$transform_variables(res1)
  z2 <- mod2$transform_variables(res2)

  z2[["prob_strain"]][, , -1] <- z2[["prob_strain"]][, 2:1, -1, drop = FALSE]
  z2[["cum_sympt_cases_non_variant_over25"]] <-
    z2[["cum_sympt_cases_over25"]] - z2[["cum_sympt_cases_non_variant_over25"]]
  z2$cum_infections_per_strain <-
    z2$cum_infections_per_strain[2:1, , drop = FALSE]
  ## This one can't easily be computed as it's not quite running
  ## incidence but over a sawtooth; the calculation relative to
  ## cum_infections_per_strain is confirmed elsewhere so here just
  ## move it out the way:
  z2[["sympt_cases_non_variant_over25_inc"]] <-
    z1[["sympt_cases_non_variant_over25_inc"]]
  for (nm in c("T_sero_neg", "R", "T_PCR_neg")) {
    z2[[nm]] <- z2[[nm]][, 2:1, , , drop = FALSE]
  }
  v5 <- c("E", "I_A", "I_P", "I_C_1", "I_C_2", "T_PCR_pre", "T_PCR_pos",
          "T_sero_pre", "T_sero_pos", "G_D", "ICU_pre_unconf", "ICU_pre_conf",
          "H_R_unconf", "H_R_conf", "H_D_unconf",
          "H_D_conf", "ICU_W_R_unconf", "ICU_W_R_conf",
          "ICU_W_D_unconf", "ICU_W_D_conf", "ICU_D_unconf",
          "ICU_D_conf", "W_R_unconf", "W_R_conf",
          "W_D_unconf", "W_D_conf")
  for (nm in v5) {
    z2[[nm]] <- z2[[nm]][, 2:1, , , , drop = FALSE]
  }

  expect_identical(z1, z2)
})


test_that("Cannot calculate Rt for multistrain without correct inputs", {
  ## Run model with 2 variants
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  expect_error(
    carehomes_Rt(steps, S[, 1, ], p),
    "Expected prob_strain input because there is more than one strain")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[-1, 1, ]),
      "Expected 'prob_strain' to have 38 rows = 19 groups x 2 strains")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, -1]),
    "Expected 'prob_strain' to have 85 columns, following 'step'")

  expect_error(
    carehomes_Rt_trajectories(steps, S, p),
    "Expected prob_strain input because there is more than one strain")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[1, , ]),
    "Expected a 3d array of 'prob_strain'")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[-1, , ]),
    "Expected 'prob_strain' to have 38 rows = 19 groups x 2 strains")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[, -1, ]),
    "Expected 2nd dim of 'prob_strain' to have length 3, following 'pars'")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[, , -1]),
    "Expected 3rd dim of 'prob_strain' to have length 85, following 'step'")

})


test_that("Can calculate Rt with an empty second variant ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ])
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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


test_that("Can calculate Rt with a second less infectious variant", {
  ## seed with 10 cases on same day as other variant

  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0.1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ])
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

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

  ## Rt should be lower (or equal) for the two variant version
  eps <- 1e-7
  expect_true(all(rt_1$Rt_all - rt_1_single_class$Rt_all <= eps))
  expect_true(all(rt_1$Rt_general - rt_1_single_class$Rt_general <= eps))
  expect_true(all(rt_all$Rt_all - rt_all_single_class$Rt_all <= eps))
  expect_true(all(rt_all$Rt_general - rt_all_single_class$Rt_general <= eps))
})


test_that("Can calculate Rt with a second more infectious variant", {
  ## Seed with 10 cases on same day as other variant
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 5),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0))


  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ])
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

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

  ## Rt should be higher (or equal) for the two variant version
  expect_true(all(rt_1$Rt_all >= rt_1_single_class$Rt_all))
  expect_true(all(rt_1$Rt_general >= rt_1_single_class$Rt_general))
  expect_true(all(rt_all$Rt_all >= rt_all_single_class$Rt_all))
  expect_true(all(rt_all$Rt_general >= rt_all_single_class$Rt_general))
})


test_that("If prob_strain is NA then Rt is NA ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  ## Remove the initial infectives so that no-one becomes infected
  info <- mod$info()
  initial <- carehomes_initial(info, 1, p)
  initial$state[info$index$I_A] <- 0

  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  ## all values of prob_strain after the first step should be NA
  expect_true(all(is.na(prob_strain[, , -1L])))

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ])
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain)

  ## all values of Rt after the first step should be NA
  expect_true(all(is.na(rt_1$Rt_all[-1L])))
  expect_true(all(is.na(rt_1$Rt_general[-1L])))
  expect_true(all(is.na(rt_1$eff_Rt_all[-1L])))
  expect_true(all(is.na(rt_1$eff_Rt_general[-1L])))
  expect_true(all(is.na(rt_all$Rt_all[-1L, ])))
  expect_true(all(is.na(rt_all$Rt_general[-1L, ])))
  expect_true(all(is.na(rt_all$eff_Rt_all[-1L, ])))
  expect_true(all(is.na(rt_all$eff_Rt_general[-1L, ])))

})


test_that("calculate Rt with both second variant and vaccination", {
  ## Seed with rate of 10 cases on same day as other variant
  ##
  ## run model with unvaccinated & vaccinated (with susceptibility halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  region <- "london"

  vacc_time <- 50
  daily_doses <- c(rep(0, vacc_time - 1), rep(Inf, vacc_time))
  n <- vaccine_priority_population(region, uptake = 1)
  vaccine_schedule <- vaccine_schedule_future(
    0, daily_doses, mean_days_between_doses = 1000, n)

  reduced_susceptibility <- 0.1 # can put anything <1 here
  transm_new_variant <- 5

  p <- carehomes_parameters(0, region,
                            waning_rate = 0,
                            strain_transmission = c(1, transm_new_variant),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            rel_susceptibility = c(1, reduced_susceptibility),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 1),
                            vaccine_progression_rate = c(0, 0),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  for (k in seq_len(np)) { # for each particle
    rt <- carehomes_Rt(steps, S[, k, ], p, prob_strain[, k, ])

    ## Impact of variant on Rt is as expected:
    ## Rt_general should increase over time because of invasion of new
    ## more transmissible variant
    ## the final value should be the initial value * the transmission advantage
    expect_approx_equal(last(rt$Rt_general),
                        rt$Rt_general[1] * transm_new_variant)

    ## Impact of vaccination on Rt is as expected:
    ## eff_Rt_general should suddenly decrease at the time at which we vaccinate
    ## everyone
    ## the reduction should be approximately by a factor reduced_susceptibility
    expect_approx_equal(rt$eff_Rt_general[vacc_time] * reduced_susceptibility,
                        rt$eff_Rt_general[vacc_time + 1],
                        rel_tol = 0.2)
  }

})


test_that("strain_rel_gamma works as expected in carehomes_parameters", {
  expect_silent(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                     strain_rel_gamma_A = 1,
                                     strain_rel_gamma_P = 1,
                                     strain_rel_gamma_C_1 = 1,
                                     strain_rel_gamma_C_2 = 1))
  expect_silent(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                     strain_rel_gamma_A = 1:2,
                                     strain_rel_gamma_P = 1:2,
                                     strain_rel_gamma_C_1 = 1:2,
                                     strain_rel_gamma_C_2 = 1:2,
                                     strain_transmission = c(1, 1)))
  expect_error(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_rel_gamma_A = c(1, 5)),
               "1 or 1")
  expect_error(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_rel_gamma_A = c(1, 5),
                                    strain_transmission = c(1, 2, 3)),
               "1 or 3")
  expect_error(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_transmission = c(1, 1),
                                    strain_rel_gamma_A = c(2, 5)),
               "must be 1")
  expect_error(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_transmission = c(1, 1),
                                    strain_rel_gamma_A = c(1, -1)),
               "non-negative")
})


test_that("carehomes_parameters_progression works as expected", {
  gammas <- c("gamma_A", "gamma_P", "gamma_C_1", "gamma_C_2")
  defaults <- c(1 / 2.88, 1 / 1.68, 1 / 2.14, 1 / 1.86)
  expect_equal(
    as.numeric(carehomes_parameters_progression(1, 1, 1, 1)[gammas]),
    defaults
  )
  expect_equal(
    as.numeric(carehomes_parameters_progression(2, 2, 2, 2)[gammas]),
    defaults * 2
  )
  expect_equal(
    matrix(unlist(carehomes_parameters_progression(1:4, 1:4, 1:4, 1:4)[gammas]),
           ncol = 4),
    vapply(defaults, function(x) x * 1:4, numeric(4))
  )
})


test_that("Relative gamma = 1 makes no difference", {
  p1 <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 1))
  p3 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 1))
  np <- 10
  mod1 <- carehomes$new(p1, 0, np, seed = 1L)
  mod2 <- carehomes$new(p2, 0, np, seed = 1L)
  mod3 <- carehomes$new(p3, 0, np, seed = 1L)
  end <- sircovid_date("2020-03-31") / p1$dt

  mod1$set_index(carehomes_index(mod1$info())$run)
  mod2$set_index(carehomes_index(mod2$info())$run)
  mod3$set_index(carehomes_index(mod3$info())$run)

  initial1 <- carehomes_initial(mod1$info(), 1, p1)
  initial2 <- carehomes_initial(mod1$info(), 1, p2)
  initial3 <- carehomes_initial(mod1$info(), 1, p3)

  res1 <- mod1$run(end)
  res2 <- mod2$run(end)
  res3 <- mod3$run(end)

  expect_equal(res2, res1)
  expect_equal(res3, res1)
})


test_that("Lower rate variant has higher Rt", {
  ## rate is .1 times ref
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = c(1, .1),
                            strain_rel_gamma_P = c(1, .1),
                            strain_rel_gamma_C_1 = c(1, .1),
                            strain_rel_gamma_C_2 = c(1, .1),
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_15_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain)

  ## rate equal to ref
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain)

  ## Rt should be higher (or equal) for the two variant version
  expect_true(all(rt_1_all$Rt_all <= rt_15_all$Rt_all))
  expect_true(all(rt_1_all$Rt_general <= rt_15_all$Rt_general))
  expect_true(all(rt_1_all$Rt_all <= rt_15_all$Rt_all))
  expect_true(all(rt_1_all$Rt_general <= rt_15$Rt_general))
})


test_that("Stuck when gamma =  0", {
  np <- 3L

  ## gammaP is 0 so IC1 is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = 1,
                            strain_rel_gamma_P = 1,
                            strain_rel_gamma_C_1 = 1,
                            strain_rel_gamma_C_2 = 1,
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0))
  p$gamma_P[] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  index_I_A <- mod$info()$index$I_A
  index_I_P <- mod$info()$index$I_P
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_2 <- mod$info()$index$I_C_2

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_true(all(unlist(y[index_I_C_1, , ]) == 0))
  expect_false(all(unlist(y[index_I_P, , ]) == 0))

  ## gammaC1 is 0 so IC2 is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = 1,
                            strain_rel_gamma_P = 1,
                            strain_rel_gamma_C_1 = 1,
                            strain_rel_gamma_C_2 = 1,
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0))
  p$gamma_C_1[] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  index_I_A <- mod$info()$index$I_A
  index_I_P <- mod$info()$index$I_P
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_2 <- mod$info()$index$I_C_2

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_true(all(unlist(y[index_I_C_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_1, , ]) == 0))

  ## gammaA is 0 & gammaC2 is 0 so R is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = 1,
                            strain_rel_gamma_P = 1,
                            strain_rel_gamma_C_1 = 1,
                            strain_rel_gamma_C_2 = 1,
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0))
  p$gamma_A[] <- 0
  p$gamma_C_2[] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  index_I_A <- mod$info()$index$I_A
  index_I_P <- mod$info()$index$I_P
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_2 <- mod$info()$index$I_C_2
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_true(all(unlist(y[index_R, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_A, , ]) == 0))
})



test_that("Stuck when gamma =  0 for second strain", {
  np <- 3L

  ## gammaP is 0 so IC1 is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = 1,
                            strain_rel_gamma_P = c(1, 0),
                            strain_rel_gamma_C_1 = 1,
                            strain_rel_gamma_C_2 = 1,
                            strain_seed_date =
                              rep(sircovid_date("2020-02-07"), 2),
                            strain_seed_rate = c(10, 0))
  p$gamma_P <- c(0, 1)

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  index_I_A <- mod$info()$index$I_A
  index_I_P <- mod$info()$index$I_P
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_2 <- mod$info()$index$I_C_2

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_true(all(unlist(y[index_I_C_1, , ]) == 0))
  expect_false(all(unlist(y[index_I_P, , ]) == 0))

  ## gammaC1 is 0 so IC2 is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = 1,
                            strain_rel_gamma_P = 1,
                            strain_rel_gamma_C_1 = 1,
                            strain_rel_gamma_C_2 = 1,
                            strain_seed_date =
                              rep(sircovid_date("2020-02-07"), 2),
                            strain_seed_rate = c(10, 0))
  p$gamma_C_1[] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  index_I_A <- mod$info()$index$I_A
  index_I_P <- mod$info()$index$I_P
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_2 <- mod$info()$index$I_C_2

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_true(all(unlist(y[index_I_C_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_1, , ]) == 0))

  ## gammaA is 0 & gammaC2 is 0 so R is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = 1,
                            strain_rel_gamma_P = 1,
                            strain_rel_gamma_C_1 = 1,
                            strain_rel_gamma_C_2 = 1,
                            strain_seed_date =
                              rep(sircovid_date("2020-02-07"), 2),
                            strain_seed_rate = c(10, 0))
  p$gamma_A[] <- 0
  p$gamma_C_2[] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

  index_I_A <- mod$info()$index$I_A
  index_I_P <- mod$info()$index$I_P
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_2 <- mod$info()$index$I_C_2
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_true(all(unlist(y[index_R, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_A, , ]) == 0))
})

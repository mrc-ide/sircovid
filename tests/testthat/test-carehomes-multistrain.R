context("carehomes (multistrain)")

test_that("carehomes_parameters_strain works as expected", {
  expect_error(
    carehomes_parameters_strain(NULL, NULL, NULL, 1),
    "At least one value required for 'strain_transmission'")
  expect_error(
    carehomes_parameters_strain(-1, NULL, NULL, 1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(c(1, -1), NULL, NULL, 1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 2), NULL, NULL, 1),
    "'strain_transmission[1]' must be 1",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 1), NULL, NULL, 1),
    "'strain_transmission[1]' must be 1",
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
    "As 'strain_seed_date' is NULL, expected 'strain_seed_value' to be NULL")
  expect_error(
    carehomes_parameters_strain(1, c(10, 20), 1, 0.25),
    "Can't use 'strain_seed_date' if only using one strain")
  expect_error(
    carehomes_parameters_strain(c(1, 1), c(10, 20, 30), 1, 0.25),
    "'strain_seed_date', if given, must be exactly two elements")
  expect_error(
    carehomes_parameters_strain(c(1, 1), c(10, 20), c(1, 1), 0.25),
    "'strain_seed_value' must be a scalar if 'strain_seed_date' is used")
  expect_error(
    carehomes_parameters_strain(c(1, 1), c(10, 20), NULL, 0.25),
    "'strain_seed_value' must be a scalar if 'strain_seed_date' is used")
})


test_that("Can seed with one-day window", {
  date <- c("2020-03-01", "2020-03-01")
  value <- 100
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), value, 1 / 4)
  expect_equal(sum(p$strain_seed_step), 100)
  expect_equal(tail(p$strain_seed_step, 6), c(0, 25, 25, 25, 25, 0))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-02"))
})


test_that("Can seed with multiple-day window", {
  date <- c("2020-03-01", "2020-03-10")
  value <- 100
  p <- carehomes_parameters_strain(c(1, 1), sircovid_date(date), value, 1 / 4)
  expect_equal(sum(p$strain_seed_step), 100 * 10)
  expect_equal(tail(p$strain_seed_step, 6), c(25, 25, 25, 25, 25, 0))
  expect_equal(sircovid_date_as_date(length(p$strain_seed_step) / 4),
               as.Date("2020-03-11"))
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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding)),
                            strain_seed_value = n_seeded_new_strain_inf)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))
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
                            strain_seed_date = c(date_seeding, date_seeding),
                            strain_seed_value = n_seeded_new_strain_inf)

  mod <- carehomes$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

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
                            strain_seed_date = c(date_seeding, date_seeding),
                            strain_seed_value = n_seeded_new_strain_inf)

  mod <- carehomes$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))
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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            waning_rate = 1 / 20,
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding)),
                            strain_seed_value = n_seeded_new_strain_inf)

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(all(y$N_tot3 - mod$transform_variables(y0)$N_tot3 == 0))
  expect_true(all(y$N_tot2 - mod$transform_variables(y0)$N_tot2 == 0))
  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot2 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot3 == 0))
})


test_that("No infection after seeding of second strain with 0 transmission", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding)),
                            strain_seed_value = n_seeded_new_strain_inf)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  ## Expect the seeded cases did not infect any other people
  expect_true(y$cum_infections_per_strain[2, 101] == n_seeded_new_strain_inf)
})


test_that("Everyone is infected when second strain transmission is large", {
  n_seeded_new_strain_inf <- 10
  date_seeding <- "2020-03-07"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1e9),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding)),
                            strain_seed_value = n_seeded_new_strain_inf)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))
  steps <- seq(0, 400, by = 4)
  date <- sircovid_date_as_date(steps / 4)
  s_date <- sircovid_date(date)
  s_date_seeding <- sircovid_date(date_seeding)
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
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

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
  y$I_asympt <- y$I_asympt[, 2:1, , , drop = FALSE]
  y$PCR_pos <- y$PCR_pos[, 2:1, , , drop = FALSE]
  y$R_pre <- y$R_pre[, 2:1, , , drop = FALSE]

  initial2_state <- unlist(y)

  mod$set_state(initial$state, initial$step)
  index <- mod$info()$index
  index_run <- c(icu = index[["I_ICU_tot"]],
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
  res1 <- dust::dust_iterate(mod, steps, index_run)

  mod2 <- carehomes$new(p, 0, np, seed = 1L, n_threads = 10)
  mod2$set_state(initial2_state, initial$step)
  mod2$set_index(index_run)
  res2 <- dust::dust_iterate(mod2, steps, index_run)

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
  y$I_asympt <- y$I_asympt[, 2:1, , , drop = FALSE]
  y$PCR_pos <- y$PCR_pos[, 2:1, , , drop = FALSE]
  y$R_pre <- y$R_pre[, 2:1, , , drop = FALSE]

  initial2_state <- unlist(y)
  mod$set_state(initial$state, initial$step)
  index <- mod$info()$index

  steps <- seq(initial$step, end, by = 1)

  res1 <- drop(dust::dust_iterate(mod, steps))

  mod2 <- carehomes$new(p, 0, np, seed = 1L)

  mod2$set_state(initial2_state, initial$step)
  res2 <- drop(dust::dust_iterate(mod2, steps))

  z1 <- mod$transform_variables(res1)
  z2 <- mod2$transform_variables(res2)

  z2$cum_infections_per_strain <-
    z2$cum_infections_per_strain[2:1, , drop = FALSE]
  for (nm in c("R_neg", "R", "PCR_neg")) {
    z2[[nm]] <- z2[[nm]][, 2:1, , , drop = FALSE]
  }
  v5 <- c("E", "I_asympt", "I_sympt", "PCR_pre", "PCR_pos", "R_pre",
          "R_pos", "I_comm_D", "I_triage_unconf", "I_triage_conf",
          "I_hosp_R_unconf", "I_hosp_R_conf", "I_hosp_D_unconf",
          "I_hosp_D_conf", "I_ICU_S_R_unconf", "I_ICU_S_R_conf",
          "I_ICU_S_D_unconf", "I_ICU_S_D_conf", "I_ICU_D_unconf",
          "I_ICU_D_conf", "R_stepdown_R_unconf", "R_stepdown_R_conf",
          "R_stepdown_D_unconf", "R_stepdown_D_conf")
  for (nm in v5) {
    z2[[nm]] <- z2[[nm]][, 2:1, , , , drop = FALSE]
  }

  expect_identical(z1, z2)
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
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- dust::dust_iterate(mod, steps, index)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

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
  y <- dust::dust_iterate(mod, steps, index)

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
                              rep(sircovid_date("2020-02-07"), 2),
                            strain_seed_value = 10)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- dust::dust_iterate(mod, steps, index)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

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
  y <- dust::dust_iterate(mod, steps, index)

  rt_1_single_class <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all_single_class <- carehomes_Rt_trajectories(steps, y, p)

  ## Effective Rt will be different because the number of infections
  ## will differ but Rt should be the same
  expect_equal(rt_1$Rt_all, rt_1_single_class$Rt_all)
  expect_equal(rt_1$Rt_general, rt_1_single_class$Rt_general)
  expect_equal(rt_all$Rt_all, rt_all_single_class$Rt_all)
  expect_equal(rt_all$Rt_general, rt_all_single_class$Rt_general)
})


test_that("Can calculate Rt with a second more infectious variant", {
  ## Seed with 10 cases on same day as other variant
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 10),
                            strain_seed_date =
                              rep(sircovid_date("2020-02-07"), 2),
                            strain_seed_value = 10)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  mod$set_index(integer(0))
  index <- mod$info()$index$S

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- dust::dust_iterate(mod, steps, index)

  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

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
  y <- dust::dust_iterate(mod, steps, index)

  rt_1_single_class <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all_single_class <- carehomes_Rt_trajectories(steps, y, p)

  ## Effective Rt will be different because the number of infections
  ## will differ but Rt should be the same with a factor of 10
  expect_equal(rt_1$Rt_all, rt_1_single_class$Rt_all * 10)
  expect_equal(rt_1$Rt_general, rt_1_single_class$Rt_general * 10)
  expect_equal(rt_all$Rt_all, rt_all_single_class$Rt_all * 10)
  expect_equal(rt_all$Rt_general, rt_all_single_class$Rt_general * 10)
})

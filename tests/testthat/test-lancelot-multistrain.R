context("lancelot (multistrain)")


test_that("lancelot_parameters_strain works as expected", {
  expect_error(
    lancelot_parameters_strain(NULL, NULL, NULL, NULL, 1),
    "At least one value required for 'strain_transmission'")
  expect_error(
    lancelot_parameters_strain(c(1, -1), NULL, NULL, NULL, 1),
    "'strain_transmission' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    lancelot_parameters_strain(c(1, 1, 1), NULL, NULL, NULL, 1),
    "Only 1 or 2")
  expect_error(
    lancelot_parameters_strain(rep(0.5, 2), NULL, NULL, NULL, 1),
    "'strain_transmission[[1]]' must be 1",
    fixed = TRUE)
  expect_error(
    lancelot_parameters_strain(rep(0.5, 1), NULL, NULL, NULL, 1),
    "'strain_transmission[[1]]' must be 1",
    fixed = TRUE)
  expect_equal(
    lancelot_parameters_strain(1, NULL, NULL, NULL, 1),
    list(n_strains = 1,
         strain_transmission = 1,
         strain_seed_step_start = 0,
         strain_seed_value = 0))
  expect_equal(
    lancelot_parameters_strain(c(1, 1), NULL, NULL, NULL, 1),
    list(n_strains = 4,
         strain_transmission = c(1, 1, 1, 1),
         strain_seed_step_start = 0,
         strain_seed_value = 0))
  expect_equal(
    lancelot_parameters_strain(c(1, 2), NULL, NULL, NULL, 1),
    list(n_strains = 4,
         strain_transmission = c(1, 2, 2, 1),
         strain_seed_step_start = 0,
         strain_seed_value = 0))
})


test_that("Prevent impossible seedings", {
  expect_error(
    lancelot_parameters_strain(c(1, 1), NULL, 1, NULL, 0.25),
    "As 'strain_seed_date' is NULL, expected 'strain_seed_size' to be NULL")
  expect_error(
    lancelot_parameters_strain(c(1, 1), NULL, NULL, 1, 0.25),
    "As 'strain_seed_date' is NULL, expected 'strain_seed_pattern' to be NULL")
  expect_error(
    lancelot_parameters_strain(1, 10, NULL, NULL, 0.25),
    "Can't use 'strain_seed_date' if only using one strain")
  expect_error(
    lancelot_parameters_strain(c(1, 1), c(10, 11), 1, 1,  1),
    "'strain_seed_date' must be a single date",
    fixed = TRUE)
  expect_error(
    lancelot_parameters_strain(c(1, 1), 10, c(1, 1), 1,  1),
    "'strain_seed_size' must be a single value",
    fixed = TRUE)
  expect_error(
    lancelot_parameters_strain(c(1, 1), 10, -1, 1,  1),
    "'strain_seed_size' must have only non-negative values",
    fixed = TRUE)
  expect_error(
    lancelot_parameters_strain(c(1, 1), 10, 1, 0,  1),
    "'strain_seed_pattern' must have only positive values",
    fixed = TRUE)
})


test_that("Can seed with one-day window", {
  date <- "2020-03-01"
  size <- 100
  pattern <- rep(1, 4)
  p <- lancelot_parameters_strain(c(1, 1), sircovid_date(date), size,
                                  pattern, 1 / 4)
  expect_equal(p$strain_seed_value, rep(25, 4))
  expect_equal(sircovid_date_as_date(p$strain_seed_step_start / 4),
               as.Date("2020-03-01"))
})


test_that("Can seed with ten-day window", {
  date <- "2020-03-01"
  size <- 100
  pattern <- rep(1, 40)
  p <- lancelot_parameters_strain(c(1, 1), sircovid_date(date), size,
                                  pattern, 1 / 4)
  expect_equal(p$strain_seed_value, rep(2.5, 40))
  expect_equal(sircovid_date_as_date(p$strain_seed_step_start / 4),
               as.Date("2020-03-01"))
})


test_that("Can seed with date which is not a multiple of step size", {
  date <- 10.2
  size <- 100
  pattern <- rep(1, 4)
  p <- lancelot_parameters_strain(c(1, 1), date, size,
                                  pattern, 1 / 4)
  expect_equal(p$strain_seed_value, c(5, 25, 25, 25, 20))
  expect_equal(sircovid_date_as_date(p$strain_seed_step_start / 4),
               sircovid_date_as_date(10))
})


test_that("Adding empty strains makes no difference", {
  p1 <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0),
                            cross_immunity = 0)
  p3 <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0.5),
                            cross_immunity = 0)
  np <- 10
  mod1 <- lancelot$new(p1, 0, np, seed = 1L)
  mod2 <- lancelot$new(p2, 0, np, seed = 1L)
  mod3 <- lancelot$new(p3, 0, np, seed = 1L)
  end <- sircovid_date("2020-03-31") / p1$dt

  mod1$set_index(lancelot_index(mod1$info())$run)
  mod2$set_index(lancelot_index(mod2$info())$run)
  mod3$set_index(lancelot_index(mod3$info())$run)

  initial1 <- lancelot_initial(mod1$info(), 1, p1)
  initial2 <- lancelot_initial(mod1$info(), 1, p2)
  initial3 <- lancelot_initial(mod1$info(), 1, p3)

  res1 <- mod1$run(end)
  res2 <- mod2$run(end)
  res3 <- mod3$run(end)

  expect_equal(res2, res1)
  expect_equal(res3, res1)
})


test_that("Seeding of second strain generates an epidemic", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  initial <- lancelot_initial(info, 1, p)
  mod$update_state(state = initial$state, step = initial$step)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  ## Did the seeded cases go on to infect other people?
  expect_true(
    all(y$cum_infections_per_strain[, 101] > n_seeded_new_strain_inf))
  ## did we count infections per strain properly?
  expect_equal(sum(y$cum_infections_per_strain[, 101]),
               y$cum_infections[, 101])

  ## Check the epidemic of the second strain starts when we expect
  steps <- seq(0, 400, by = 4)
  date <- sircovid_date_as_date(steps / 4)
  s_date <- sircovid_date(date)
  s_date_seeding <- sircovid_date(date_seeding)

  ## Expect to see cases recorded after seeding
  ## No cases before seeding
  expect_true(all(y$E[, 2:4, , , s_date < s_date_seeding + 1] == 0))
  ## No cases on seeding day other than in 4th age group
  expect_true(all(y$E[-4, 2:4, , , s_date == s_date_seeding + 1] == 0))
  ## Some cases on seeding day in 4th age group
  expect_true(y$E[4, 2, 1, , s_date == s_date_seeding + 1] > 0)
  ## No seeding into strains 3 and 4
  expect_true(all(y$E[4, 3:4, 1, , s_date == s_date_seeding + 1] == 0))
  ## Some cases on all days after seeding day

  ## It's not guaranteed that *all* will be greater than zero, but most will be
  ## Tolerance decreased to 0.123 to account for time to get to 3 and 4
  expect_true(
    mean(colSums(y$E[, 2:4, 1, , s_date >= s_date_seeding]) > 0) > 0.123)
})


test_that("Second more virulent strain takes over", {
  np <- 10
  n_seeded_new_strain_inf <- 10
  start_date <- sircovid_date("2020-02-07")
  date_seeding <- start_date # seed both strains on same day
  p <- lancelot_parameters(start_date, "england",
                           strain_transmission = c(1, 10),
                           strain_seed_date = date_seeding,
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))

  mod <- lancelot$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## cumulative infections with 2nd strain larger than with 1st strain
  ## (average over 10 runs)
  expect_true(mean(y$cum_infections_per_strain[c(1, 3), , 101]) <
                mean(y$cum_infections_per_strain[c(2, 4), , 101]))
})


test_that("Second less virulent strain does not take over", {
  np <- 10
  n_seeded_new_strain_inf <- 10
  start_date <- sircovid_date("2020-02-07")
  date_seeding <- start_date # seed both strains on same day
  p <- lancelot_parameters(start_date, "england",
                           strain_transmission = c(1, 0.1),
                           strain_seed_date = date_seeding,
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))

  mod <- lancelot$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  ## Cumulative infections with 2nd strain smaller than with 1st strain
  ## (average over 10 runs)
  expect_true(mean(y$cum_infections_per_strain[c(1, 3), , 101]) >
                mean(y$cum_infections_per_strain[c(2, 4), , 101]))

})


test_that("N_tots stay constant with second strain and no waning immunity - no
          superinfection", {
  ## Default for waning_rate is 0
  set.seed(1)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))

  mod <- lancelot$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(y$N_tot_sero_1 -
                    mod$transform_variables(y0)$N_tot_sero_1 == 0))
  expect_true(all(y$N_tot_sero_2 -
                    mod$transform_variables(y0)$N_tot_sero_2 == 0))
  expect_true(all(y$N_tot_PCR - mod$transform_variables(y0)$N_tot_PCR == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot_sero_1 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot_sero_2 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot_PCR == 0))
})


test_that("N_tot is constant with second strain and waning immunity, while
          sero N_tots are non-decreasing - superinfection", {
  ## Default for waning_rate is 0
  set.seed(1)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           waning_rate = 1 / 20,
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(diff(y$N_tot_sero_1) >= 0))
  expect_true(all(diff(y$N_tot_sero_2) >= 0))
  expect_true(all(diff(y$N_tot_PCR) >= 0))
  expect_true(all(colSums(y$N_tot) <= y$N_tot_sero_1))
  expect_true(all(colSums(y$N_tot) <= y$N_tot_sero_2))
  expect_true(all(colSums(y$N_tot) <= y$N_tot_PCR))
})

test_that("N_tot is constant with second strain and waning immunity, while
          sero and PCR N_tots are non-decreasing - no superinfection", {
            ## Default for waning_rate is 0
  set.seed(1)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           waning_rate = 1 / 20,
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(diff(y$N_tot_sero_1) >= 0))
  expect_true(all(diff(y$N_tot_sero_2) >= 0))
  expect_true(all(diff(y$N_tot_PCR) >= 0))
  expect_true(all(colSums(y$N_tot) <= y$N_tot_sero_1))
  expect_true(all(colSums(y$N_tot) <= y$N_tot_sero_2))
  expect_true(all(colSums(y$N_tot) <= y$N_tot_PCR))
})


test_that("No-one in strains 3 or 4 if waning_rate is 1e6", {
  set.seed(2L)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           waning_rate = 1e6,
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, 1, seed = 2L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$R[, 3:4, , ] == 0))
})


test_that("No-one in strains 3 or 4 if no super infection", {
  set.seed(2L)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           waning_rate = 0.1,
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))

  mod <- lancelot$new(p, 0, 1, seed = 2L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$R[, 3:4, , ] == 0))
})


test_that("prob_strain sums to 1", {
  np <- 3L
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"

  # no waning immunity
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  mod <- lancelot$new(p, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_prob_strain <- mod$info()$index$prob_strain

  expect_equal(initial$state[index_prob_strain], c(1, 0))

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)
  prob_strain <- y[index_prob_strain, , ]

  expect_equal(prob_strain[1, , ], 1 - prob_strain[2, , ])

  # with waning immunity
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           waning_rate = 1 / 20,
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  mod <- lancelot$new(p, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_prob_strain <- mod$info()$index$prob_strain
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)
  prob_strain <- y[index_prob_strain, , ]

  expect_equal(prob_strain[1, , ], 1 - prob_strain[2, , ])
  expect_equal(dim(prob_strain), c(2, 3, 123))
})


test_that("No infection after seeding of second strain with 0 transmission", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 0),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## can't compare to fixed number so generate a feasible range
  pois_range <- rpois(1000, n_seeded_new_strain_inf)

  ## Expect the seeded cases did not infect any other people
  expect_true(y$cum_infections_per_strain[2, 101] <= max(pois_range))
  expect_true(y$cum_infections_per_strain[2, 101] >= min(pois_range))
  expect_true(y$cum_infections_per_strain[4, 101] <= max(pois_range))
  expect_true(y$cum_infections_per_strain[4, 101] >= 0)
})


test_that("Everyone is infected when second strain transmission is large", {
  n_seeded_new_strain_inf <- 10
  date_seeding <- "2020-03-07"
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1e9),
                           strain_seed_date = sircovid_date(date_seeding),
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  ## set gamma_E to Inf so that seeded individuals move through each E stage
  ## in one step
  p$gamma_E_step <- Inf

  mod <- lancelot$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  steps <- seq(0, 400, by = 4)
  date <- sircovid_date_as_date(steps / 4)
  s_date <- sircovid_date(date)
  s_date_seeding <- sircovid_date(date_seeding[[1]])
  ## No cases before seeding
  expect_true(all(y$E[, 2:4, , , s_date < s_date_seeding] == 0))

  ## The +2 is because we need seeded individuals to get out of the first and
  ## second E compartments before they can go on to infect others
  expect_true(all(y$S[, 1, s_date > (s_date_seeding + 2)] == 0))
})


test_that("No infection with either strain with perfect vaccine", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england",
                           strain_transmission = c(1, 1),
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           waning_rate = 0.1,
                           cross_immunity = 0)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  state <- lancelot_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  index_E <- array(info$index$E, info$dim$E)
  state[index_E[4, 2, 1, 1]] <- 10 # seed infections with second strain

  mod$update_state(state = state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## Noone moves into unvaccinated
  ## except in the group where infections because of waning immunity
  expect_true(all(y$S[-4, 1, ] == 0))

  ## Noone changes compartment within the vaccinated individuals
  expect_true(all(y$S[, 2, ] == y$S[, 2, 1]))

  ## Noone gets infected with either strain
  expect_true(all(y$cum_infections_per_strain == 0))
  expect_true(all(y$cum_infections == 0))
})


test_that("different strains are equivalent", {
  set.seed(1)
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           cross_immunity = 0)
  np <- 10
  mod <- lancelot$new(p, 0, np, seed = 1L, n_threads = 10)
  end <- sircovid_date("2020-07-31") / p$dt
  end <- sircovid_date("2020-03-31") / p$dt

  initial <- lancelot_initial(mod$info(), 1, p)
  y <- mod$transform_variables(initial$state)
  i <- c(2, 1, 4, 3)
  y$I_A <- y$I_A[, i, , , drop = FALSE]
  y$T_PCR_pos <- y$T_PCR_pos[, i, , , drop = FALSE]
  y$T_sero_pre_1 <- y$T_sero_pre_1[, i, , , drop = FALSE]
  y$T_sero_pre_2 <- y$T_sero_pre_2[, i, , , drop = FALSE]

  initial2_state <- unlist(y)

  mod$update_state(state = initial$state, step = initial$step)
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

  mod2 <- lancelot$new(p, 0, np, seed = 1L, n_threads = 10)
  mod2$update_state(state = initial2_state, step = initial$step)
  mod2$set_index(index_run)
  res2 <- mod2$simulate(steps)

  expect_equal(res1, res2)
})


test_that("Swapping strains gives identical results with different index", {
  p <- lancelot_parameters(0, "england",
                           initial_seed_size = 0,
                           strain_transmission = c(1, 1),
                           cross_immunity = 0)

  np <- 1
  mod <- lancelot$new(p, 0, np, seed = 1L)
  end <- sircovid_date("2020-05-1") / p$dt
  initial <- lancelot_initial(mod$info(), 1, p)
  y <- mod$transform_variables(initial$state)
  ## force some seeding
  y$I_A[4, 1, 1, 1] <- 10
  initial_state <- unlist(y)

  i <- c(2, 1, 4, 3)
  y$I_A <- y$I_A[, i, , , drop = FALSE]
  y$I_weighted <- y$I_weighted[, i, , drop = FALSE]
  y$prob_strain <- y$prob_strain[c(2, 1)]

  initial2_state <- unlist(y)
  mod$update_state(state = initial_state, step = initial$step)
  index <- mod$info()$index

  steps <- seq(initial$step, end, by = 1)

  res1 <- drop(mod$simulate(steps))

  mod2 <- lancelot$new(p, 0, np, seed = 1L)

  mod2$update_state(state = initial2_state, step = initial$step)
  res2 <- drop(mod2$simulate(steps))

  z1 <- mod$transform_variables(res1)
  z2 <- mod2$transform_variables(res2)

  z2[["prob_strain"]] <- z2[["prob_strain"]][2:1, , drop = FALSE]
  z2[["cum_sympt_cases_non_variant"]] <-
    z2[["cum_sympt_cases"]] - z2[["cum_sympt_cases_non_variant"]]
  z2[["cum_sympt_cases_non_variant_over25"]] <-
    z2[["cum_sympt_cases_over25"]] - z2[["cum_sympt_cases_non_variant_over25"]]
  z2$cum_infections_per_strain <-
    z2$cum_infections_per_strain[i, , drop = FALSE]
  ## This one can't easily be computed as it's not quite running
  ## incidence but over a sawtooth; the calculation relative to
  ## cum_infections_per_strain is confirmed elsewhere so here just
  ## move it out the way:
  z2[["sympt_cases_non_variant_inc"]] <-
    z1[["sympt_cases_non_variant_inc"]]
  z2[["sympt_cases_non_variant_over25_inc"]] <-
    z1[["sympt_cases_non_variant_over25_inc"]]
  for (nm in c("T_sero_neg_1", "T_sero_neg_2", "R", "T_PCR_neg",
               "I_weighted")) {
    z2[[nm]] <- z2[[nm]][, i, , , drop = FALSE]
  }
  v5 <- c("E", "I_A", "I_P", "I_C_1", "I_C_2", "T_PCR_pre", "T_PCR_pos",
          "T_sero_pre_1", "T_sero_pre_2", "T_sero_pos_1", "T_sero_pos_2",
          "G_D", "ICU_pre_unconf", "ICU_pre_conf",
          "H_R_unconf", "H_R_conf", "H_D_unconf",
          "H_D_conf", "ICU_W_R_unconf", "ICU_W_R_conf",
          "ICU_W_D_unconf", "ICU_W_D_conf", "ICU_D_unconf",
          "ICU_D_conf", "W_R_unconf", "W_R_conf",
          "W_D_unconf", "W_D_conf")
  for (nm in v5) {
    z2[[nm]] <- z2[[nm]][, i, , , , drop = FALSE]
  }

  expect_identical(z1, z2)
})


test_that("Cannot calculate Rt for multistrain without correct inputs", {
  ## Run model with 2 variants
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  ## Check lancelot_Rt R
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p),
    "Expected prob_strain and R input because")
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[-1, 1, ]),
    "Expected 'R' to have 76 rows")
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, -1]),
    "Expected 'R' to have 123 columns")

  ## Check lancelot_Rt prob_strain
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ]),
    "Expected prob_strain and R input because")
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[-1, 1, ], R = R[, 1, ]),
    "Expected a 2 strains")
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ][1, , drop = FALSE],
                R = R[, 1, ]),
    "Expected 'prob_strain' to have 2 rows")
  expect_error(
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, -1], R = R[, 1, ]),
    "Expected 'prob_strain' to have 123 columns, following 'step'")

  ## Check lancelot_Rt_trajectories R
  expect_error(
    lancelot_Rt_trajectories(steps, S, p),
    "Expected prob_strain and R input because")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain[1, , ], R = R[1, , ]),
    "Expected a 3d array of 'R'")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R[-1, , ]),
    "Expected 'R' to have 76 rows")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R[, -1, ]),
    "Expected 2nd and 3rd")

  ## Check lancelot_Rt_trajectories prob_strain
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, R = R),
    "Expected prob_strain and R input because")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain[1, , ], R = R),
    "Expected a 3d array of 'prob_strain'")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain[-1, , ], R = R),
    "Expected a 3d array of 'prob_strain'")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain[, -1, ], R = R),
    "Expected 2nd dim of 'prob_strain' to have length 3, following 'pars'")
  expect_error(
    lancelot_Rt_trajectories(steps, S, p, prob_strain[, , -1], R = R),
    "Expected 3rd dim of 'prob_strain' to have length 123, following 'step'")
})

test_that("Cannot calculate IFR_t for multistrain without correct inputs", {
  ## Run model with 2 variants
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  I <- y[index_I, , ]
  R <- y[index_R, , ]

  ## Check lancelot_ifr_t R
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p),
    "Expected R input because")
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p, R = R[-1, 1, ]),
    "Expected 'R' to have 76 rows")
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p, R = R[, 1, -1]),
    "Expected 'R' to have 123 columns")

  p1 <- p
  p1$n_strains <- 3
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p1, R = R[, 1, ]),
    "Multstrain IFR currently only works if n_strains is 4")

  ## Check lancelot_ifr_t_trajectories R
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I, p),
    "Expected R input because")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I, p, R = R[1, , ]),
    "Expected a 3d array of 'R'")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I, p, R = R[-1, , ]),
    "Expected 'R' to have 76 rows")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I, p, R = R[, -1, ]),
    "Expected 'S' and 'R' to have same length of 2nd dim")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I, p, R = R[, , -1]),
    "Expected 3rd dim of 'R' to have length 123, given 'step'")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, 1, , drop = FALSE],
                                I[, 1, , drop = FALSE], list(p), R = R),
    "Expected 2nd dim of 'R' to have length 1, given 'pars'")

})

## Tests for basic object properties, not analytical results from calculations
test_that("wtmean_Rt works as expected", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           cross_immunity = 0)
  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ])
  expect_equal(dim(rt$eff_Rt_all), c(123, 2))
  expect_equal(class(rt), c("multi_strain", "Rt"))

  rt_traj <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R)
  expect_equal(dim(rt_traj$eff_Rt_all), c(123, 2, 3))
  expect_equal(class(rt_traj), c("multi_strain", "Rt_trajectories", "Rt"))

  rt_strain_weighted <-
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                weight_Rt = TRUE)
  expect_equal(length(rt_strain_weighted$eff_Rt_all), 123)
  expect_equal(class(rt_strain_weighted), c("single_strain", "Rt"))

  rt_traj_strain_weighted <-
    lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                             weight_Rt = TRUE)
  expect_equal(dim(rt_traj_strain_weighted$eff_Rt_all), c(123, 3))
  expect_equal(class(rt_traj_strain_weighted),
               c("single_strain", "Rt_trajectories", "Rt"))

  nms <- names(rt)

  avg_rt <- wtmean_Rt(rt, prob_strain[, 1, ])
  avg_rt_traj <- wtmean_Rt(rt_traj, prob_strain)

  ## here the 1st strain has weight 1 all along
  ## (except step 1 --> to investigate)
  ## so expect the average R to be the same as R for strain 1
  expect_equal(rt$eff_Rt_general[, 1], avg_rt$eff_Rt_general)
  expect_equal(rt$eff_Rt_all[, 1], avg_rt$eff_Rt_all)
  expect_equal(rt$Rt_general[, 1], avg_rt$Rt_general)
  expect_equal(rt$Rt_all[, 1], avg_rt$Rt_all)

  expect_equal(rt_traj$eff_Rt_general[, 1, ],
               avg_rt_traj$eff_Rt_general[, ])
  expect_equal(rt_traj$eff_Rt_all[, 1, ], avg_rt_traj$eff_Rt_all[, ])
  expect_equal(rt_traj$Rt_general[, 1, ], avg_rt_traj$Rt_general[, ])
  expect_equal(rt_traj$Rt_all[, 1, ], avg_rt_traj$Rt_all[, ])

  ## the average should be the same if calculated inside the Rt calculation
  ## functions or post hoc
  expect_equal(rt_strain_weighted$eff_Rt_general,
               avg_rt$eff_Rt_general)
  expect_equal(rt_strain_weighted$eff_Rt_all,
               avg_rt$eff_Rt_all)
  expect_equal(rt_strain_weighted$Rt_general,
               avg_rt$Rt_general)
  expect_equal(rt_strain_weighted$Rt_all,
               avg_rt$Rt_all)

  expect_equal(names(avg_rt), nms)
  expect_equal(names(avg_rt_traj), nms)

  expect_error(wtmean_Rt(1L), "must inherit")
  expect_error(wtmean_Rt(structure(list(Rt_all = matrix(1)), class = "Rt"),
                         prob_strain),
               "Expect elements of Rt to have dimensions")

  ## check single particle case
  S <- mcstate::array_flatten(S, 2:3)[, 1, drop = FALSE]
  R <- mcstate::array_flatten(R, 2:3)[, 1, drop = FALSE]
  prob_strain <- mcstate::array_flatten(prob_strain, 2:3)[, 1, drop = FALSE]
  rt_weight_F <- lancelot_Rt(1, S, p, prob_strain, R = R, weight_Rt = FALSE)
  rt_weight_T <- lancelot_Rt(1, S, p, prob_strain, R = R, weight_Rt = TRUE)
  expect_equal(rt_weight_F$eff_Rt_all[[1]], rt_weight_T$eff_Rt_all)
  expect_equal(rt_weight_F$eff_Rt_general[[1]],
               rt_weight_T$eff_Rt_general)
  expect_equal(rt_weight_F$Rt_all[[1]], rt_weight_T$Rt_all)
  expect_equal(rt_weight_F$Rt_general[[1]], rt_weight_T$Rt_general)
})


test_that("Can calculate Rt with an empty second variant ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  expect_equal(rt_1$step, rt_1_single_class$step)
  expect_equal(rt_1$date, rt_1_single_class$date)
  expect_equal(rt_1$beta, rt_1_single_class$beta)
  expect_equal(rt_1$eff_Rt_all, rt_1_single_class$eff_Rt_all)
  expect_equal(rt_1$eff_Rt_general,
               rt_1_single_class$eff_Rt_general)
  expect_equal(rt_1$Rt_all, rt_1_single_class$Rt_all)
  expect_equal(rt_1$Rt_general, rt_1_single_class$Rt_general)

  expect_equal(rt_all$step, rt_all_single_class$step)
  expect_equal(rt_all$date, rt_all_single_class$date)
  expect_equal(rt_all$beta, rt_all_single_class$beta)
  expect_equal(rt_all$eff_Rt_all,
               rt_all_single_class$eff_Rt_all)
  expect_equal(rt_all$eff_Rt_general,
               rt_all_single_class$eff_Rt_general)
  expect_equal(rt_all$Rt_all, rt_all_single_class$Rt_all)
  expect_equal(rt_all$Rt_general,
               rt_all_single_class$Rt_general)
})


test_that("Can calculate Rt with strain_transmission (1, 0) ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 0),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = FALSE)
  expect_vector_equal(rt_1$eff_Rt_all[, 2], 0)
  expect_vector_equal(rt_1$eff_Rt_general[, 2], 0)
  expect_vector_equal(rt_1$Rt_all[, 2], 0)
  expect_vector_equal(rt_1$Rt_general[, 2], 0)
})


test_that("Can calculate Rt with a second less infectious variant", {
  ## seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 0.1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 2L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 2L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be lower (or equal) for the two variant version
  tol <- 1e-5
  expect_vector_lte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_lte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_lte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_lte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("Can calculate Rt with a second more infectious variant", {
  ## Seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 5),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)


  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be higher (or equal) for the two variant version
  tol <- 1e-5
  expect_vector_gte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_gte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("Can calculate Rt with a second less lethal variant", {
  ## Seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_severity = c(1, 0),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)


  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be higher (or equal) for the two variant version
  ## because less letal --> more people recover and they have longer
  ## duration of infection
  tol <- 1e-5
  expect_vector_gte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_gte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("Can calculate Rt with a second variant with longer I_A", {
  ## Seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = c(1, 0.1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)


  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be higher (or equal) for the two variant version
  ## because longer duration of infection (for asymptomatics)
  ## added the 0.001 as seems to be rounding error issues
  tol <- 1e-3
  expect_vector_gte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_gte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("Can calculate Rt with a second variant with longer I_P", {
  ## Seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_P = c(1, 0.1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)


  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be higher (or equal) for the two variant version
  ## because longer duration of infection (for presymptomatics)
  ## added the 0.001 as seems to be rounding error issues
  tol <- 1e-3
  expect_vector_gte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_gte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("Can calculate Rt with a second variant with longer I_C_1", {
  ## Seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_C_1 = c(1, 0.1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)


  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- lancelot_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- lancelot_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be higher (or equal) for the two variant version
  ## because longer duration of infection (for I_C_1)
  ## added the 0.001 as seems to be rounding error issues
  tol <- 1e-3
  expect_vector_gte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_gte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("If prob_strain is NA then Rt is NA only in same steps", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1))

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  info <- mod$info()
  initial <- lancelot_initial(info, 1, p)

  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  ## set NA for prob_strain in steps 60-70
  na_steps <- 60:70
  prob_strain[, , na_steps] <- NA

  expect_true(!any(is.na(prob_strain[, , -na_steps])))
  expect_true(all(is.na(prob_strain[, , na_steps])))

  ## test unweighted

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ])
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R)

  ## not all NA as weight_Rt = FALSE
  expect_vector_equal(lengths(rt_1[1:3]), 123)
  expect_true(!any(is.na(simplify2array(rt_1[4:7])[na_steps, , ])))
  expect_true(!any(is.na(simplify2array(rt_1[4:7])[-na_steps, , ])))

  expect_equal(dim(simplify2array(rt_all[1:3])), c(123, 3, 3))
  expect_true(!any(is.na(simplify2array(rt_all[4:7])[na_steps, , , ])))
  expect_true(!any(is.na(simplify2array(rt_all[4:7])[-na_steps, , , ])))


  ## test weighted

  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## all values of Rt in 60:70 to be NA and others not to be
  expect_vector_equal(lengths(rt_1[1:3]), 123)
  expect_true(all(is.na(simplify2array(rt_1[4:7])[na_steps, ])))
  expect_true(!any(is.na(simplify2array(rt_1[4:7])[-na_steps, ])))

  expect_equal(dim(simplify2array(rt_all[1:3])), c(123, 3, 3))
  expect_true(all(is.na(simplify2array(rt_all[4:7])[na_steps, , ])))
  expect_true(!any(is.na(simplify2array(rt_all[4:7])[-na_steps, , ])))
})


test_that("Can compute Rt if all prob_strain is NA", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1))

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  info <- mod$info()
  initial <- lancelot_initial(info, 1, p)

  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  ## set all NA for prob_strain
  prob_strain[] <- NA

  ## test unweighted
  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ])
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R)

  ## No NA as weight_Rt = FALSE
  expect_vector_equal(vlapply(rt_1, function(x) any(is.na(x))), FALSE)
  expect_vector_equal(vlapply(rt_all, function(x) any(is.na(x))), FALSE)

  ## test weighted
  rt_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                      weight_Rt = TRUE)
  rt_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                     weight_Rt = TRUE)

  ## All NA as weight_Rt = TRUE
  expect_vector_equal(vlapply(rt_1[1:3], function(x) any(is.na(x))), FALSE)
  expect_vector_equal(vlapply(rt_1[4:7], function(x) all(is.na(x))), TRUE)
  expect_vector_equal(vlapply(rt_all[1:3], function(x) any(is.na(x))), FALSE)
  expect_vector_equal(vlapply(rt_all[4:7], function(x) all(is.na(x))), TRUE)
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

  reduced_susceptibility <- 1 # can put anything <1 here
  transm_new_variant <- 5

  p <- lancelot_parameters(0, region,
                           waning_rate = 0.1,
                           strain_transmission = c(1, transm_new_variant),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           rel_susceptibility = rep(1, 3),
                           vaccine_progression_rate = numeric(3),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L,
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 2L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain
  index_I <- mod$info()$index$I_weighted

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(2)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  I <- y[index_I, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  for (k in seq_len(np)) { # for each particle

    rt <- lancelot_Rt(steps, S[, k, ], p, prob_strain[, k, ], R = R[, k, ],
                      weight_Rt = TRUE)

    ## Impact of variant on Rt is as expected:
    ## Rt_general should increase over time because of invasion of new
    ## more transmissible variant
    ## the final value should be the initial value * the transmission advantage
    expect_approx_equal(last(rt$Rt_general),
                        rt$Rt_general[1] * transm_new_variant)

    ## Impact of vaccination on Rt is as expected:
    ## eff_Rt_general should suddenly decrease at the time at which we
    ## vaccinate everyone the reduction should be approximately by a factor
    ## reduced_susceptibility
    expect_approx_equal(rt$eff_Rt_general[vacc_time] * reduced_susceptibility,
                        rt$eff_Rt_general[vacc_time + 1],
                        rel_tol = 0.2)

    ## Not sure how we want to test results here yet, but just make sure it runs
    ifr_t <- lancelot_ifr_t(steps, S[, k, ], I[, k, ], p, R = R[, k, ])
  }

})

test_that("strain_rel_severity works as expected in lancelot_parameters", {
  strain_rel_severity <- c(1, 0.5)
  rel_p_death <- c(1, 0.6, 0.7)
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_severity = strain_rel_severity,
                           rel_p_death = rel_p_death)
  # check strains are mirrored
  expect_equal(p$rel_p_ICU_D[, 1:2, ], p$rel_p_ICU_D[, 4:3, ])
  expect_equal(p$rel_p_ICU_D[, 2, ],
               p$rel_p_ICU_D[, 1, ] * strain_rel_severity[2])
  expect_equal(p$rel_p_ICU_D[, , 2],
               p$rel_p_ICU_D[, , 1] * rel_p_death[2])
  expect_equal(p$rel_p_ICU_D[, , 3],
               p$rel_p_ICU_D[, , 1] * rel_p_death[3])

})

test_that("strain_rel_gamma works as expected in lancelot_parameters", {
  expect_silent(lancelot_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_rel_gamma_A = 1,
                                    strain_rel_gamma_P = 1,
                                    strain_rel_gamma_C_1 = 1,
                                    strain_rel_gamma_C_2 = 1))
  expect_silent(lancelot_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_rel_gamma_A = 1:2,
                                    strain_rel_gamma_P = 1:2,
                                    strain_rel_gamma_C_1 = 1:2,
                                    strain_rel_gamma_C_2 = 1:2,
                                    strain_transmission = c(1, 1)))
  expect_error(lancelot_parameters(sircovid_date("2020-02-07"), "england",
                                   strain_rel_gamma_A = c(1, 5)),
               "1 or 1")
  expect_error(lancelot_parameters(sircovid_date("2020-02-07"), "england",
                                   strain_rel_gamma_A = c(1, 5),
                                   strain_transmission = c(1, 2, 3)),
               "1 or 2")
  expect_error(lancelot_parameters(sircovid_date("2020-02-07"), "england",
                                   strain_transmission = c(1, 1),
                                   strain_rel_gamma_A = c(2, 5)),
               "must be 1")
  expect_error(lancelot_parameters(sircovid_date("2020-02-07"), "england",
                                   strain_transmission = c(1, 1),
                                   strain_rel_gamma_A = c(1, -1)),
               "non-negative")
})


test_that("Relative gamma = 1 makes no difference", {
  p1 <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            cross_immunity = 0)
  p3 <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            cross_immunity = 0)
  np <- 10
  mod1 <- lancelot$new(p1, 0, np, seed = 1L)
  mod2 <- lancelot$new(p2, 0, np, seed = 1L)
  mod3 <- lancelot$new(p3, 0, np, seed = 1L)
  end <- sircovid_date("2020-03-31") / p1$dt

  mod1$set_index(lancelot_index(mod1$info())$run)
  mod2$set_index(lancelot_index(mod2$info())$run)
  mod3$set_index(lancelot_index(mod3$info())$run)

  initial1 <- lancelot_initial(mod1$info(), 1, p1)
  initial2 <- lancelot_initial(mod1$info(), 1, p2)
  initial3 <- lancelot_initial(mod1$info(), 1, p3)

  res1 <- mod1$run(end)
  res2 <- mod2$run(end)
  res3 <- mod3$run(end)

  expect_equal(res2, res1)
  expect_equal(res3, res1)
})


test_that("Lower rate variant has higher Rt", {
  ## rate is .1 times ref
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = c(1, .1),
                           strain_rel_gamma_P = c(1, .1),
                           strain_rel_gamma_C_1 = c(1, .1),
                           strain_rel_gamma_C_2 = c(1, .1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_15_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                        weight_Rt = TRUE)

  ## rate equal to ref
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_1_all <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                       weight_Rt = TRUE)

  ## Rt should be higher (or equal) for the two variant version
  expect_vector_lte(rt_1_all$Rt_all, rt_15_all$Rt_all)
  expect_vector_lte(rt_1_all$Rt_general, rt_15_all$Rt_general)
})


test_that("Stuck when gamma =  0", {
  np <- 3L

  ## gammaP is 0 so IC1 is 0
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = 1,
                           strain_rel_gamma_P = 1,
                           strain_rel_gamma_C_1 = 1,
                           strain_rel_gamma_C_2 = 1,
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  p$gamma_P_step <- 0

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)

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
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = 1,
                           strain_rel_gamma_P = 1,
                           strain_rel_gamma_C_1 = 1,
                           strain_rel_gamma_C_2 = 1,
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  p$gamma_C_1_step <- 0

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)

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
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = 1,
                           strain_rel_gamma_P = 1,
                           strain_rel_gamma_C_1 = 1,
                           strain_rel_gamma_C_2 = 1,
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4))
  p$gamma_A_step <- 0
  p$gamma_C_2_step <- 0

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)

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

  ## gammaP is 0 so IC1 is 0 for second strain
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = c(1, 1),
                           strain_rel_gamma_P = c(1, 0),
                           strain_rel_gamma_C_1 = c(1, 1),
                           strain_rel_gamma_C_2 = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)

  strain_1 <- c(1:19, 58:76)
  strain_2 <- c(20:38, 39:57)

  index_I_A <- mod$info()$index$I_A
  index_I_A_strain_2 <- index_I_A[strain_2]
  index_I_A_strain_1 <- index_I_A[strain_1]
  index_R <- mod$info()$index$R
  index_R_strain_2 <- index_R[strain_2]
  index_R_strain_1 <- index_R[strain_1]
  index_I_P <- mod$info()$index$I_P
  index_I_P_strain_1 <- index_I_P[strain_1]
  index_I_P_strain_2 <- index_I_P[strain_2]
  index_I_C_1 <- mod$info()$index$I_C_1
  index_I_C_1_strain_1 <- index_I_C_1[strain_1]
  index_I_C_1_strain_2 <- index_I_C_1[strain_2]
  index_I_C_2 <- mod$info()$index$I_C_2
  index_I_C_2_strain_1 <- index_I_C_2[strain_1]
  index_I_C_2_strain_2 <- index_I_C_2[strain_2]

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  expect_false(all(unlist(y[index_I_C_1_strain_1, , ]) == 0))
  expect_true(all(unlist(y[index_I_C_1_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_P_strain_2, , ]) == 0))

  ## gammaC1 is 0 so IC2 is 0 for second strain
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = c(1, 1),
                           strain_rel_gamma_P = c(1, 1),
                           strain_rel_gamma_C_1 = c(1, 0),
                           strain_rel_gamma_C_2 = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, np, seed = 1L)
  mod$update_state(state = initial$state, step = initial$step)
  set.seed(1)
  y <- mod$simulate(steps)
  expect_false(all(unlist(y[index_I_C_2_strain_1, , ]) == 0))
  expect_true(all(unlist(y[index_I_C_2_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_1_strain_2, , ]) == 0))


  ## gammaA is 0 & gammaC2 is 0 so R is 0
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_gamma_A = c(1, 0),
                           strain_rel_gamma_P = c(1, 1),
                           strain_rel_gamma_C_1 = c(1, 1),
                           strain_rel_gamma_C_2 = c(1, 0),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, np, seed = 1L)
  mod$update_state(state = initial$state, step = initial$step)
  set.seed(1)
  y <- mod$simulate(steps)
  expect_false(all(unlist(y[index_R_strain_1, , ]) == 0))
  expect_true(all(unlist(y[index_R_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_2_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_2_strain_2, , ]) == 0))
})


test_that("No one is hospitalised, no-one recovers in edge case 2 - multi", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           waning_rate = 1 / 20,
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 1,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  p$p_C_step[, ] <- 1
  p$p_H_step[, ] <- 1
  p$p_G_D_step[, ] <- 1

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  ## Move initial infectives to 2nd stage sympt
  y0 <- lancelot_initial(info, 1, p)$state
  y0[info$index$I_C_2] <- y0[info$index$I_A]
  y0[info$index$I_A] <- 0

  mod$update_state(state = y0)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_C_2 > 0))
  expect_true(all(y$H_R_unconf == 0))
  expect_true(all(y$H_R_conf == 0))
  expect_true(all(y$H_D_unconf == 0))
  expect_true(all(y$H_D_conf == 0))
  expect_true(all(y$ICU_W_R_unconf == 0))
  expect_true(all(y$ICU_W_R_conf == 0))
  expect_true(all(y$ICU_W_D_unconf == 0))
  expect_true(all(y$ICU_W_D_conf == 0))
  expect_true(all(y$ICU_D_unconf == 0))
  expect_true(all(y$ICU_D_conf == 0))
  expect_true(all(y$ICU_pre_unconf == 0))
  expect_true(all(y$ICU_pre_conf == 0))
  expect_true(all(y$W_R_unconf == 0))
  expect_true(all(y$W_R_conf == 0))
  expect_true(all(y$W_D_unconf == 0))
  expect_true(all(y$W_D_conf == 0))
  expect_true(all(y$R == 0))
  expect_true(all(y$D_hosp == 0))
})


test_that("G_D empty when p_G_D = 0", {
  np <- 3L
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 1,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  p$p_G_D_step[, ] <- 0

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)

  expect_true(all(y[mod$info()$index$G_D, , ] == 0))
})


test_that("G_D strain 2 empty when p_G_D = c(1, 0)", {
  np <- 3L
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_severity = c(1, 0),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  # Strain 1 and 4 (2 -> 1)
  expect_false(all(y$G_D[, 1, , , , ] == 0))
  expect_false(all(y$G_D[, 4, , , , ] == 0))
  # Strain 2 and 3 (1 -> 2)
  expect_true(all(y$G_D[, 2, , , , ] == 0))
  expect_true(all(y$G_D[, 3, , , , ] == 0))
})


test_that("Can't move from S to E3/4", {
  np <- 3L
  start_date <- sircovid_date("2020-02-07")
  p <- lancelot_parameters(start_date, "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))
  ## we are seeding
  expect_true(any(y$E[, 1, , , , start_date + 2] > 0))
  expect_true(any(y$E[, 2, , , , start_date + 2] > 0))
  expect_true(all(y$E[, 3:4, , , , start_date + 2] == 0))
})


test_that("Nobody in R2-R4 when strain_transmission = c(1, 0)", {
  np <- 3L
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 0),
                           beta_value = 1e100,
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_true(all(y$S[, , , 123] == 0))
  ## there will be the seeded individuals in the 15-19 age band
  expect_true(all(y$R[-4, 2:4, , , 123] == 0))
  expect_false(all(y$R[, 1, , , 123] == 0))
})


test_that("Can only move to S from R3 and R4 to S", {
  np <- 3L
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 0,
                           strain_transmission = c(1, 1),
                           waning_rate = 1 / 5,
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  ## Prevent anyone leaving S
  p$rel_susceptibility[] <- 0

  mod <- lancelot$new(p, 0, np, seed = 1L)
  info <- mod$info()

  initial <- lancelot_initial(info, np, p)
  y0 <- initial$state
  ## Empty R1 and R2
  y0[info$index$R][1:38] <- 0
  ## Fill R3 and R4
  y0[info$index$R][39:76] <- 1e3

  mod$update_state(state = y0)

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  diff_S <- y$S[-4, , , 2] - y$S[-4, , , 1]
  diff_R <- y$R[-4, 3:4, , , 1] - y$R[-4, 3:4, , , 2]
  diff_R <- apply(diff_R, c(1, 3), sum)
  expect_equal(diff_R, diff_S)
  expect_true(all(y$E[-4, , , , , ] == 0))
})


test_that("Everyone in R3 and R4 when no waning and transmission high", {
  np <- 1L
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 30,
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 30,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  p$strain_transmission[] <- 1e8
  ## set p_C to 0 so that individuals move to R quickly
  p$p_C_step[, ] <- 0

  mod <- lancelot$new(p, 0, np, seed = 1L)
  info <- mod$info()

  initial <- lancelot_initial(mod$info(), np, p)
  y0 <- initial$state

  mod$update_state(state = y0, step = initial$step)
  end <- sircovid_date("2021-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_true(all(y$R[, 1:2, , 450] == 0))
  expect_false(all(y$R[, 3:4, , 450] == 0))
})


test_that("cross_immunity parameter errors when expected", {
  expect_error(
    lancelot_parameters(sircovid_date("2020-02-07"), "england",
                        cross_immunity = 2), "in [0, 1]", fixed = TRUE
  )
  expect_error(
    lancelot_parameters(sircovid_date("2020-02-07"), "england",
                        cross_immunity = -2), "in [0, 1]", fixed = TRUE
  )
  expect_error(
    lancelot_parameters(sircovid_date("2020-02-07"), "england",
                        cross_immunity = c(1, 1)), "Invalid length"
  )
})

test_that("complete cross_immunity means no Strain 3/4 infections", {
  np <- 1L

  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 1)
  mod <- lancelot$new(p, 0, np)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_true(all(y$cum_infections_per_strain[3:4, ] == 0))
})

test_that("some cross-immunity means less Strain 3 or 4 infections than none
           and > 0", {
 skip_on_windows_gha() # see #356
 np <- 1L

 ## no cross-immnunity
 p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                          strain_transmission = c(1, 1),
                          strain_seed_date = sircovid_date("2020-02-07"),
                          strain_seed_size = 10,
                          strain_seed_pattern = rep(1, 4),
                          cross_immunity = 0)
 mod <- lancelot$new(p, 0, np)
 initial <- lancelot_initial(mod$info(), np, p)
 mod$update_state(state = initial$state, step = initial$step)
 end <- sircovid_date("2020-05-01") / p$dt
 steps <- seq(initial$step, end, by = 1 / p$dt)
 set.seed(1)
 y <- mod$transform_variables(
   drop(mod$simulate(steps)))
 infect_no_cross <- y$cum_infections_per_strain[3:4, 123]

 ## some cross-immunity
 p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                          strain_transmission = c(1, 1),
                          strain_seed_date = sircovid_date("2020-02-07"),
                          strain_seed_size = 10,
                          strain_seed_pattern = rep(1, 4),
                          cross_immunity = 0.5)
 mod <- lancelot$new(p, 0, np)
 initial <- lancelot_initial(mod$info(), np, p)
 mod$update_state(state = initial$state, step = initial$step)
 set.seed(1)
 y <- mod$transform_variables(
   drop(mod$simulate(steps)))
 infect_some_cross <- y$cum_infections_per_strain[3:4, 123]

 expect_vector_lte(infect_some_cross, infect_no_cross)
 expect_true(all(infect_some_cross > 0))
 expect_true(all(infect_no_cross > 0))
})


test_that("cross-immunity can be separated by strain", {
  seed <- 1
  np <- 1L
  set.seed(seed)

  ## complete immunity from Strain 1 means Strain 3 empty (1 -> 2)
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 100,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = c(1, 0)
  )
  mod <- lancelot$new(p, 0, np, seed = seed)
  initial <- lancelot_initial(mod$info(), np, p)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  mod$update_state(state = initial$state, step = initial$step)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_equal(y$cum_infections_per_strain[3, 123], 0)
  expect_gt(y$cum_infections_per_strain[4, 123], 0)

  ## complete immunity from Strain 2 means Strain 4 empty (2 -> 1)
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 100,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = c(0, 1))
  mod <- lancelot$new(p, 0, np, seed = seed)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_equal(y$cum_infections_per_strain[4, 123], 0)
  expect_gt(y$cum_infections_per_strain[3, 123], 0)
})


test_that("Can calculate ifr_t with an empty second variant ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           rel_susceptibility = c(1, 1, 1),
                           cross_immunity = 0)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  I <- y[index_I, , ]
  R <- y[index_R, , ]

  ifr_t_1 <- lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p, R = R[, 1, ])
  ifr_t_all <- lancelot_ifr_t_trajectories(steps, S, I, p, R = R)

  ## Run model with one strain only
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           rel_susceptibility = c(1, 1, 1))

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_I))
  y <- mod$simulate(steps)

  S <- y[1:57, , ]
  I <- y[58:114, , ]

  ifr_t_1_single_class <- lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p)
  ifr_t_all_single_class <- lancelot_ifr_t_trajectories(steps, S, I, p)

  expect_equal(ifr_t_1$step, ifr_t_1_single_class$step)
  expect_equal(ifr_t_1$date, ifr_t_1_single_class$date)
  expect_equal(ifr_t_1$IFR_t_all, ifr_t_1_single_class$IFR_t_all)
  expect_equal(ifr_t_1$IFR_t_general, ifr_t_1_single_class$IFR_t_general)
  expect_equal(ifr_t_1$IHR_t_all, ifr_t_1_single_class$IHR_t_all)
  expect_equal(ifr_t_1$IHR_t_general, ifr_t_1_single_class$IHR_t_general)
  expect_equal(ifr_t_1$IFR_t_all_no_vacc,
               ifr_t_1_single_class$IFR_t_all_no_vacc)
  expect_equal(ifr_t_1$IFR_t_general_no_vacc,
               ifr_t_1_single_class$IFR_t_general_no_vacc)
  expect_equal(ifr_t_1$IHR_t_all_no_vacc,
               ifr_t_1_single_class$IHR_t_all_no_vacc)
  expect_equal(ifr_t_1$IHR_t_general_no_vacc,
               ifr_t_1_single_class$IHR_t_general_no_vacc)
  expect_equal(ifr_t_1$ALOS, ifr_t_1_single_class$ALOS)
  expect_equal(ifr_t_1$ALOS_no_vacc, ifr_t_1_single_class$ALOS_no_vacc)

  expect_equal(ifr_t_all$step, ifr_t_all_single_class$step)
  expect_equal(ifr_t_all$date, ifr_t_all_single_class$date)
  expect_equal(ifr_t_all$IFR_t_all, ifr_t_all_single_class$IFR_t_all)
  expect_equal(ifr_t_all$IFR_t_general, ifr_t_all_single_class$IFR_t_general)
  expect_equal(ifr_t_all$IHR_t_all, ifr_t_all_single_class$IHR_t_all)
  expect_equal(ifr_t_all$IHR_t_general, ifr_t_all_single_class$IHR_t_general)
  expect_equal(ifr_t_all$IFR_t_all_no_vacc,
               ifr_t_all_single_class$IFR_t_all_no_vacc)
  expect_equal(ifr_t_all$IFR_t_general_no_vacc,
               ifr_t_all_single_class$IFR_t_general_no_vacc)
  expect_equal(ifr_t_all$IHR_t_all_no_vacc,
               ifr_t_all_single_class$IHR_t_all_no_vacc)
  expect_equal(ifr_t_all$IHR_t_general_no_vacc,
               ifr_t_all_single_class$IHR_t_general_no_vacc)
  expect_equal(ifr_t_all$ALOS, ifr_t_all_single_class$ALOS)
  expect_equal(ifr_t_all$ALOS_no_vacc, ifr_t_all_single_class$ALOS_no_vacc)
})


test_that("Can calculate ifr_t with a second less lethal variant", {
  ## Seed with 10 cases on same day as other variant
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_rel_severity = c(1, 0.5),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 1)


  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  I <- y[index_I, , ]
  R <- y[index_R, , ]

  ifr_t_1 <- lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p, R = R[, 1, ])
  ifr_t_all <- lancelot_ifr_t_trajectories(steps, S, I, p, R = R)

  ## move everyone into first strain
  I[1:19, , ] <- I[1:19, , ] + I[20:38, , ] + I[39:57, , ] + I[58:76, , ]
  I[20:76, , ] <- 0
  R[1:19, , ] <- R[1:19, , ] + R[20:38, , ] + R[39:57, , ] + R[58:76, , ]
  R[20:76, , ] <- 0

  ifr_t_1_empty_strain2 <-
    lancelot_ifr_t(steps, S[, 1, ], I[, 1, ], p, R = R[, 1, ])
  ifr_t_all_empty_strain2 <-
    lancelot_ifr_t_trajectories(steps, S, I, p, R = R)

  ## Rt should be lower (or equal) for the two variant version
  ## because less lethal
  tol <- 1e-5
  expect_vector_lte(ifr_t_1$IFR_t_all, ifr_t_1_empty_strain2$IFR_t_all,
                    tol = tol)
  expect_vector_lte(ifr_t_1$IFR_t_general, ifr_t_1_empty_strain2$IFR_t_general,
                    tol = tol)
  expect_vector_lte(ifr_t_all$IFR_t_all, ifr_t_all_empty_strain2$IFR_t_all,
                    tol = tol)
  expect_vector_lte(ifr_t_all$IFR_t_general,
                    ifr_t_all_empty_strain2$IFR_t_general, tol = tol)
})


test_that("can inflate the number of strains after running with 1", {
  ## one strain
  p1 <- lancelot_parameters(sircovid_date("2020-02-07"), "england")
  ## two strains
  p2 <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1))

  np <- 3L
  mod1 <- lancelot$new(p1, 0, np, seed = 1L)
  initial <- lancelot_initial(mod1$info(), 10, p1)
  mod1$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p1$dt
  steps <- seq(initial$step, end, by = 1 / p1$dt)
  y1 <- mod1$run(end)
  info1 <- mod1$info()

  mod2 <- lancelot$new(p2, 0, 1)
  info2 <- mod2$info()
  y2 <- inflate_state_strains(y1, info1, info2)

  expect_equal(dim(y2), c(info2$len, 3))
  expect_equal(sum(y2), sum(y1))

  expect_equal(sort(y2[, 1], decreasing = TRUE)[seq_len(info1$len)],
               sort(y1[, 1], decreasing = TRUE))
  expect_equal(sort(y2[, 1], decreasing = TRUE)[-seq_len(info1$len)],
               rep(0, info2$len - info1$len))

  ## examples of the different types of conversions:
  ## rank 1, length 2:
  z1 <- mod1$transform_variables(y1)
  z2 <- mod2$transform_variables(y2)
  expect_equal(z2$prob_strain, matrix(c(1, 0), 2, 3))
  expect_equal(z2$cum_infections_per_strain,
               rbind(z1$cum_infections_per_strain, 0, 0, 0))
  expect_equal(z2$T_PCR_pos[, 1, , , , drop = FALSE], z1$T_PCR_pos)
  expect_equal(z2$T_PCR_pos[, 2:4, , , , drop = FALSE],
               array(0, c(19, 3, 2, 1, 3)))
  expect_equal(z2$R[, 1, , , drop = FALSE], z1$R)
  expect_equal(z2$R[, 2:4, , , drop = FALSE],
               array(0, c(19, 3, 1, 3)))

  expect_equal(z2$time, z1$time)

  expect_error(inflate_state_strains(1), "Expected a matrix")
  expect_error(inflate_state_strains(matrix(1), list(len = 2), 2),
               "Expected a matrix with 2 rows")
  expect_error(
    inflate_state_strains(matrix(1, 2, 2),
                          list(len = 2, index = list(a = 1)),
                          list(len = 2, index = list(b = 1))),
    "Can't inflate state")
})


test_that("Rt lower with perfect cross immunity", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 0.1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 1)

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt_cross_1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ],
                            R = R[, 1, ], weight_Rt = TRUE)

  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 0.1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)

  rt_cross_0 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ],
                            R = R[, 1, ], weight_Rt = TRUE)

  ## Rt should be equal for strain 1
  tol <- 1e-5
  expect_vector_lte(rt_cross_1$Rt_all, rt_cross_0$Rt_all, tol = tol)
  expect_vector_lte(rt_cross_1$Rt_general, rt_cross_0$Rt_general, tol = tol)
  expect_vector_lt(rt_cross_1$eff_Rt_all, rt_cross_0$eff_Rt_all, tol = tol)
  expect_vector_lt(rt_cross_1$eff_Rt_general, rt_cross_0$eff_Rt_general,
                   tol = tol)
})


test_that("Can calculate Rt with asymmetric cross immunity", {
  ## Run with perfect cross immunity for both strains first
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 10,
                           initial_seed_pattern = rep(1, 4),
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = c(1, 1))

  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$update_state(state = initial$state, step = initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  rt <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ],
                    R = R[, 1, ], weight_Rt = FALSE)

  ## Now calculate where strain 1 gives perfect immunity to strain 2,
  ## and strain 2 gives no immunity to strain 1
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 10,
                           initial_seed_pattern = rep(1, 4),
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = c(1, 0))

  rt1 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ],
                     R = R[, 1, ], weight_Rt = FALSE)

  ## Rt should be unaffected
  tol <- 1e-5
  expect_equal(rt1$Rt_all, rt$Rt_all)
  expect_equal(rt1$Rt_general, rt$Rt_general)
  ## eff_Rt should be unaffected for strain 2, increased for strain 1
  expect_equal(rt1$eff_Rt_all[, 2], rt$eff_Rt_all[, 2], tol = tol)
  expect_equal(rt1$eff_Rt_general[, 2], rt$eff_Rt_general[, 2], tol = tol)
  expect_vector_gt(rt1$eff_Rt_all[, 1], rt$eff_Rt_all[, 1], tol = tol)
  expect_vector_gt(rt1$eff_Rt_general[, 1], rt$eff_Rt_general[, 1], tol = tol)


  ## Now calculate the other way round; where strain 2 gives perfect immunity
  ## to strain 1, and strain 1 gives no immunity to strain 2
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           initial_seed_size = 10,
                           initial_seed_pattern = rep(1, 4),
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = c(0, 1))

  rt2 <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ],
                     R = R[, 1, ], weight_Rt = FALSE)

  ## Rt should be unaffected
  tol <- 1e-5
  expect_equal(rt2$Rt_all, rt$Rt_all)
  expect_equal(rt2$Rt_general, rt$Rt_general)
  ## eff_Rt should be unaffected for strain 1, increased for strain 2
  expect_equal(rt2$eff_Rt_all[, 1], rt$eff_Rt_all[, 1], tol = tol)
  expect_equal(rt2$eff_Rt_general[, 1], rt$eff_Rt_general[, 1], tol = tol)
  expect_vector_gt(rt2$eff_Rt_all[, 2], rt$eff_Rt_all[, 2], tol = tol)
  expect_vector_gt(rt2$eff_Rt_general[, 2], rt$eff_Rt_general[, 2], tol = tol)

})


test_that("Can interpolate multistrain Rt", {
  dat <- reference_data_lancelot_mcmc()
  rt <- local({
    p <- lapply(seq_len(nrow(dat$pars)), function(i)
      dat$predict$transform(dat$pars[i, ]))
    i <- grep("S_", rownames(dat$trajectories$state))
    S <- dat$trajectories$state[i, , ]
    lancelot_Rt_trajectories(dat$trajectories$step, S, p)
  })

  future <- list(
    "2020-04-15" = future_Rt(1.5, "2020-03-10"),
    "2020-05-01" = future_Rt(0.5, "2020-03-10"),
    "2020-05-15" = future_Rt(2, "2020-03-10"))


  res <- future_relative_beta(future, rt$date[, 1], rt$Rt_general)
  baseline <- add_future_betas(dat, rt, future)

  events <- sircovid_simulate_events("2020-03-30", "2020-06-01", NULL)
  p <- lapply(seq_len(nrow(baseline$pars)), function(i)
    baseline$predict$transform(baseline$pars[i, ]))
  ans <- sircovid_simulate(lancelot, baseline$state, p, events,
                           index = dat$predict$index)

  ## Work out our critical dates so that we can start interpolation:
  step <- attr(ans, "step")
  S <- ans[grep("S_", rownames(ans)), , ]

  crit_dates <- sircovid_date(names(future))

  set.seed(1)
  rt_cmp <- lancelot_Rt_trajectories(step, S, p,
                                     initial_step_from_parameters = FALSE)

  ## Only interpolate if "every" is given:
  set.seed(1)
  expect_identical(
    lancelot_Rt_trajectories(step, S, p,
                             initial_step_from_parameters = FALSE,
                             interpolate_min = 3),
    rt_cmp)


  ## Then compute the Rt values with interpolation
  rt_int_2 <- lancelot_Rt_trajectories(step, S, p,
                                       initial_step_from_parameters = FALSE,
                                       interpolate_every = 2,
                                       interpolate_min = 3,
                                       interpolate_critical_dates = crit_dates)
  rt_int_7 <- lancelot_Rt_trajectories(step, S, p,
                                       initial_step_from_parameters = FALSE,
                                       interpolate_every = 7,
                                       interpolate_min = 3,
                                       interpolate_critical_dates = crit_dates)
  rt_int_14 <- lancelot_Rt_trajectories(step, S, p,
                                        initial_step_from_parameters = FALSE,
                                        interpolate_every = 14,
                                        interpolate_min = 1,
                                        interpolate_critical_dates =
                                          crit_dates)
  ## check the error is small
  tol <- 0.05
  # for interpolation every 2 days
  expect_vector_equal(rt_cmp$eff_Rt_all, rt_int_2$eff_Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$eff_Rt_general, rt_int_2$eff_Rt_general,
                      tol = tol)
  expect_vector_equal(rt_cmp$Rt_all, rt_int_2$Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$Rt_general, rt_int_2$Rt_general, tol = tol)
  # for interpolation every 7 days
  expect_vector_equal(rt_cmp$eff_Rt_all, rt_int_7$eff_Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$eff_Rt_general,
                      rt_int_7$eff_Rt_general, tol = tol)
  expect_vector_equal(rt_cmp$Rt_all, rt_int_7$Rt_all, tol = tol)
  expect_vector_equal(rt_cmp$Rt_general, rt_int_7$Rt_general, tol = tol)
  # have to increase tolerance dramatically for every 14 days
  tol2 <- 0.5
  expect_vector_equal(rt_cmp$eff_Rt_all, rt_int_14$eff_Rt_all, tol = tol2)
  expect_vector_equal(rt_cmp$eff_Rt_general,
                      rt_int_14$eff_Rt_general, tol = tol2)
  expect_vector_equal(rt_cmp$Rt_all, rt_int_14$Rt_all, tol = tol2)
  expect_vector_equal(rt_cmp$Rt_general, rt_int_14$Rt_general, tol = tol2)
})



test_that("wtmean_Rt works as expected with interpolation", {
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           cross_immunity = 0)
  np <- 3L
  mod <- lancelot$new(p, 0, np, seed = 1L)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]
  interpolate_every <- 7
  interpolate_min <- 3

  rt <- lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                    interpolate_every = interpolate_every,
                    interpolate_min = interpolate_min)
  expect_equal(dim(rt$eff_Rt_all), c(123, 2))
  expect_equal(class(rt), c("multi_strain", "Rt"))

  rt_traj <- lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                      interpolate_every = interpolate_every,
                                      interpolate_min = interpolate_min)
  expect_equal(dim(rt_traj$eff_Rt_all), c(123, 2, 3))
  expect_equal(class(rt_traj), c("multi_strain", "Rt_trajectories", "Rt"))

  rt_strain_weighted <-
    lancelot_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                weight_Rt = TRUE, interpolate_every = interpolate_every,
                interpolate_min = interpolate_min)
  expect_equal(length(rt_strain_weighted$eff_Rt_all), 123)
  expect_equal(class(rt_strain_weighted), c("single_strain", "Rt"))

  rt_traj_strain_weighted <-
    lancelot_Rt_trajectories(steps, S, p, prob_strain, R = R,
                             weight_Rt = TRUE,
                             interpolate_every = interpolate_every,
                             interpolate_min = interpolate_min)
  expect_equal(dim(rt_traj_strain_weighted$eff_Rt_all), c(123, 3))
  expect_equal(class(rt_traj_strain_weighted),
               c("single_strain", "Rt_trajectories", "Rt"))

  nms <- names(rt)

  avg_rt <- wtmean_Rt(rt, prob_strain[, 1, ])
  avg_rt_traj <- wtmean_Rt(rt_traj, prob_strain)

  ## here the 1st strain has weight 1 all along
  ## (except step 1 --> to investigate)
  ## so expect the average R to be the same as R for strain 1
  expect_equal(rt$eff_Rt_general[, 1], avg_rt$eff_Rt_general)
  expect_equal(rt$eff_Rt_all[, 1], avg_rt$eff_Rt_all)
  expect_equal(rt$Rt_general[, 1], avg_rt$Rt_general)
  expect_equal(rt$Rt_all[, 1], avg_rt$Rt_all)

  expect_equal(rt_traj$eff_Rt_general[, 1, ],
               avg_rt_traj$eff_Rt_general[, ])
  expect_equal(rt_traj$eff_Rt_all[, 1, ], avg_rt_traj$eff_Rt_all[, ])
  expect_equal(rt_traj$Rt_general[, 1, ], avg_rt_traj$Rt_general[, ])
  expect_equal(rt_traj$Rt_all[, 1, ], avg_rt_traj$Rt_all[, ])

  ## the average should be the same if calculated inside the Rt calculation
  ## functions or post hoc
  expect_equal(rt_strain_weighted$eff_Rt_general,
               avg_rt$eff_Rt_general)
  expect_equal(rt_strain_weighted$eff_Rt_all,
               avg_rt$eff_Rt_all)
  expect_equal(rt_strain_weighted$Rt_general,
               avg_rt$Rt_general)
  expect_equal(rt_strain_weighted$Rt_all,
               avg_rt$Rt_all)

  expect_equal(names(avg_rt), nms)
  expect_equal(names(avg_rt_traj), nms)

  ## check single particle case
  S <- mcstate::array_flatten(S, 2:3)[, 1, drop = FALSE]
  R <- mcstate::array_flatten(R, 2:3)[, 1, drop = FALSE]
  prob_strain <- mcstate::array_flatten(prob_strain, 2:3)[, 1, drop = FALSE]
  rt_weight_F <- lancelot_Rt(1, S, p, prob_strain, R = R, weight_Rt = FALSE)
  rt_weight_T <- lancelot_Rt(1, S, p, prob_strain, R = R, weight_Rt = TRUE)
  expect_equal(rt_weight_F$eff_Rt_all[[1]], rt_weight_T$eff_Rt_all)
  expect_equal(rt_weight_F$eff_Rt_general[[1]],
               rt_weight_T$eff_Rt_general)
  expect_equal(rt_weight_F$Rt_all[[1]], rt_weight_T$Rt_all)
  expect_equal(rt_weight_F$Rt_general[[1]], rt_weight_T$Rt_general)
})


test_that("Can rotate strains", {
  np <- 6
  n_seeded_new_strain_inf <- 10
  start_date <- sircovid_date("2020-02-07")
  date_seeding <- start_date # seed both strains on same day
  p <- lancelot_parameters(start_date, "england",
                           strain_transmission = c(1, 10),
                           strain_seed_date = date_seeding,
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))

  mod <- lancelot$new(p, 0, np, seed = 1L)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)$state
  mod$update_state(state = lancelot_initial(info, 1, p)$state)
  invisible(mod$run(200))

  state1 <- mod$state()
  state2 <- rotate_strains(state1, mod$info())

  expect_equal(dim(state2), dim(state1))
  expect_equal(colSums(state2), colSums(state1))

  y1 <- mod$transform_variables(state1)
  y2 <- mod$transform_variables(state2)

  ## Next, check one of each rank of transformed variable to make sure
  ## that all are appropriately transformed.

  ## This one is doable analytically :)
  expect_equal(y2$prob_strain, matrix(1:0, 2, np))

  expect_true(all(y2$I_C_1[, 2:4, , , ] == 0))
  expect_equal(y2$I_C_1[, 1, , , ], apply(y1$I_C_1, c(1, 5), sum))

  expect_true(all(y2$R[, 2:4, , ] == 0))
  expect_equal(y2$R[, 1, , ], apply(y1$R, c(1, 4), sum))

  expect_true(all(y2$cum_infections_per_strain[2:4, ] == 0))
  expect_equal(y2$cum_infections_per_strain[1, ],
               colSums(y1$cum_infections_per_strain))

  expect_error(rotate_strains(state1[, 1], mod$info()),
               "Expected a matrix or array for 'state'")
  expect_error(rotate_strains(state1[-1, ], mod$info()),
               "Expected a matrix with [0-9]+ rows for 'state'")

  ## Structured case; show that it does not matter in what order we
  ## structure the particle dimensions.
  state3 <- rotate_strains(mcstate::array_reshape(state1, 2, 2:3), mod$info())
  expect_equal(dim(state3), c(nrow(state1), 2, 3))
  expect_equal(state3[, 1, 1], state2[, 1])
  expect_equal(state3, mcstate::array_reshape(state2, 2, 2:3))
})


test_that("Can rotate strains with cross-infection", {
  np <- 2L

  ## no cross-immnunity
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = sircovid_date("2020-02-07"),
                           strain_seed_size = 10,
                           strain_seed_pattern = rep(1, 4),
                           cross_immunity = 0)
  mod <- lancelot$new(p, 0, np)
  initial <- lancelot_initial(mod$info(), np, p)
  mod$update_state(state = initial$state, step = initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  set.seed(1)
  invisible(mod$run(end))

  state1 <- mod$state()
  state2 <- rotate_strains(state1, mod$info())

  expect_equal(dim(state2), dim(state1))
  expect_equal(colSums(state2), colSums(state1))

  y1 <- mod$transform_variables(state1)
  y2 <- mod$transform_variables(state2)

  ## Next, check one of each rank of transformed variable to make sure
  ## that all are appropriately transformed; see above.
  expect_equal(y2$prob_strain, matrix(1:0, 2, np))
  expect_true(all(y2$I_C_1[, 2:4, , , ] == 0))
  expect_equal(y2$I_C_1[, 1, , , ], apply(y1$I_C_1, c(1, 5), sum))
  expect_true(all(y2$R[, 2:4, , ] == 0))
  expect_equal(y2$R[, 1, , ], apply(y1$R, c(1, 4), sum))
  expect_true(all(y2$cum_infections_per_strain[2:4, ] == 0))
  expect_equal(y2$cum_infections_per_strain[1, ],
               colSums(y1$cum_infections_per_strain))
})


test_that("Can rotate strains with vaccination", {
  p <- lancelot_parameters(0, "england",
                           strain_transmission = c(1, 1),
                           rel_susceptibility = c(0.5, 0.5),
                           rel_p_sympt = c(0.5, 0.5),
                           rel_p_hosp_if_sympt = c(0.5, 0.5),
                           waning_rate = 0.1,
                           cross_immunity = 0)
  np <- 3
  mod <- lancelot$new(p, 0, np, seed = 42L)
  info <- mod$info()

  state <- lancelot_initial(info, 1, p)$state
  index_E <- array(info$index$E, info$dim$E)
  state[index_E[4, , 1, 1]] <- 100

  mod$update_state(state = state)
  invisible(mod$run(200))

  state1 <- mod$state()
  state2 <- rotate_strains(state1, mod$info())

  expect_equal(dim(state2), dim(state1))
  expect_equal(colSums(state2), colSums(state1))

  y1 <- mod$transform_variables(state1)
  y2 <- mod$transform_variables(state2)

  ## Next, check one of each rank of transformed variable to make sure
  ## that all are appropriately transformed; see above.
  expect_equal(y2$prob_strain, matrix(1:0, 2, np))
  expect_true(all(y2$I_C_1[, 2:4, , , ] == 0))
  expect_equal(y2$I_C_1[, 1, , , ], apply(y1$I_C_1, c(1, 4, 5), sum))
  expect_true(all(y2$R[, 2:4, , ] == 0))
  expect_equal(y2$R[, 1, , ], apply(y1$R, c(1, 3, 4), sum))
  expect_true(all(y2$cum_infections_per_strain[2:4, ] == 0))
  expect_equal(y2$cum_infections_per_strain[1, ],
               colSums(y1$cum_infections_per_strain))
})


test_that("rotate strain uses correct variables", {
  path <- sircovid_file("odin/lancelot.R")
  json <- odin::odin_parse_(path, options = odin.dust::odin_dust_options())
  ir <- odin::odin_ir_deserialise(json)

  check1 <- function(v) {
    any(vlapply(unlist(c(v$dimnames$length, v$dimnames$dim)), function(x)
      any(c("n_strains", "n_real_strains") %in% ir$equations[[x]]$rhs$value)))
  }

  vars <- ir$data$elements[names(ir$data$variable$contents)]
  n_strain_dim <- vlapply(vars, check1)

  ## Descriptive names to get good test failures:
  variables_missing_from_rotate <-
    setdiff(names(which(n_strain_dim)), rotate_strain_compartments)
  expect_length(variables_missing_from_rotate, 0)

  variables_rotated_but_no_strain <-
    setdiff(rotate_strain_compartments, names(which(n_strain_dim)))
  expect_length(variables_rotated_but_no_strain, 0)
})

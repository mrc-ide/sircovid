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
    carehomes_parameters_strain(c(1, 1, 1), NULL, NULL, 1),
    "Only 1 or 2")
  expect_error(
    carehomes_parameters_strain(rep(0.5, 2), NULL, NULL, 1),
    "'strain_transmission[[1]]' must be 1",
    fixed = TRUE)
  expect_error(
    carehomes_parameters_strain(rep(0.5, 1), NULL, NULL, 1),
    "'strain_transmission[[1]]' must be 1",
    fixed = TRUE)
  expect_equal(
    carehomes_parameters_strain(1, NULL, NULL, 1),
    list(n_strains = 1,
         strain_transmission = 1,
         strain_seed_step = 0))
  expect_equal(
    carehomes_parameters_strain(c(1, 1), NULL, NULL, 1),
    list(n_strains = 4,
         strain_transmission = c(1, 1, 1, 1),
         strain_seed_step = 0))
  expect_equal(
    carehomes_parameters_strain(c(1, 2), NULL, NULL, 1),
    list(n_strains = 4,
         strain_transmission = c(1, 2, 2, 1),
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
                             strain_transmission = c(1, 0),
                             cross_immunity = 0)
  p3 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 0.5),
                             cross_immunity = 0)
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
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
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

  ## No cases before seeding
  expect_true(all(y$E[, 2:4, , , s_date < s_date_seeding] == 0))
  ## No cases on seeding day other than in 4th age group
  expect_true(all(y$E[-4, 2:4, , , s_date == s_date_seeding] == 0))
  ## Some cases on seeding day in 4th age group
  expect_true(y$E[4, 2, 1, , s_date == s_date_seeding] > 0)
  ## No seeding into strains 3 and 4
  expect_true(all(y$E[4, 3:4, 1, , s_date == s_date_seeding] == 0))
  ## Some cases on all days after seeding day

  ## It's not guaranteed that *all* will be greater than zero, but most will be
  ## Tolerance decreased to 0.85 to account for time to get to 3 and 4
  expect_true(
    mean(colSums(y$E[, 2:4, 1, , s_date >= s_date_seeding]) > 0) > 0.85)
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
  expect_true(mean(y$cum_infections_per_strain[c(1, 3), , 101]) <
                mean(y$cum_infections_per_strain[c(2, 4), , 101]))
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
  expect_true(mean(y$cum_infections_per_strain[c(1, 3), , 101]) >
                mean(y$cum_infections_per_strain[c(2, 4), , 101]))

})


test_that("N_tots stay constant with second strain and no waning immunity - no
          superinfection", {
  ## Default for waning_rate is 0
  set.seed(1)
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
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            waning_rate = 1 / 20,
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
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
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            waning_rate = 1e6,
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, 1, seed = 2L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$R[, 3:4, , ] == 0))
})


test_that("No-one in strains 3 or 4 if no super infection", {
  set.seed(2L)
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            waning_rate = 0.1,
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0))

  mod <- carehomes$new(p, 0, 1, seed = 2L)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$R[, 3:4, , ] == 0))
})


test_that("prob_strain sums to 1", {
  np <- 3L
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  date_seeding_end <- "2020-03-08"

  # no waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)
  mod <- carehomes$new(p, 0, np, seed = 1L)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  index_prob_strain <- mod$info()$index$prob_strain

  expect_equal(initial$state[index_prob_strain], c(1, 0))

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)
  prob_strain <- y[index_prob_strain, , ]

  expect_equal(prob_strain[1, , ], 1 - prob_strain[2, , ])

  # with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            waning_rate = 1 / 20,
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)
  mod <- carehomes$new(p, 0, np, seed = 1L)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  index_prob_strain <- mod$info()$index$prob_strain
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)
  prob_strain <- y[index_prob_strain, , ]

  expect_equal(prob_strain[1, , ], 1 - prob_strain[2, , ])
  expect_equal(dim(prob_strain), c(2, 3, 85))
})


test_that("No infection after seeding of second strain with 0 transmission", {
  n_seeded_new_strain_inf <- 100
  date_seeding <- "2020-03-07"
  date_seeding_end <- "2020-03-08"
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0),
                            strain_seed_date =
                              sircovid_date(c(date_seeding, date_seeding_end)),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)

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
  expect_true(y$cum_infections_per_strain[4, 101] <= max(pois_range))
  expect_true(y$cum_infections_per_strain[4, 101] >= 0)
})


test_that("Everyone is infected when second strain transmission is large", {
  n_seeded_new_strain_inf <- 10
  date_seeding <- c("2020-03-07", "2020-03-08")
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1e9),
                            strain_seed_date =
                              sircovid_date(date_seeding),
                            strain_seed_rate = c(n_seeded_new_strain_inf, 0),
                            cross_immunity = 0)

  ## set gamma_E to Inf so that seeded individuals move through each E stage
  ## in one step
  p$gamma_E_step <- Inf

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
  expect_true(all(y$E[, 2:4, , , s_date < s_date_seeding] == 0))

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
                            waning_rate = 0.1,
                            cross_immunity = 0)
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
  expect_true(all(y$cum_infections == 0))
})


test_that("different strains are equivalent", {
  set.seed(1)
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            cross_immunity = 0)
  np <- 10
  mod <- carehomes$new(p, 0, np, seed = 1L, n_threads = 10)
  end <- sircovid_date("2020-07-31") / p$dt
  end <- sircovid_date("2020-03-31") / p$dt

  initial <- carehomes_initial(mod$info(), 1, p)
  y <- mod$transform_variables(initial$state)
  i <- c(2, 1, 4, 3)
  y$I_A <- y$I_A[, i, , , drop = FALSE]
  y$T_PCR_pos <- y$T_PCR_pos[, i, , , drop = FALSE]
  y$T_sero_pre_1 <- y$T_sero_pre_1[, i, , , drop = FALSE]
  y$T_sero_pre_2 <- y$T_sero_pre_2[, i, , , drop = FALSE]

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
                            strain_transmission = c(1, 1),
                            cross_immunity = 0)

  np <- 1
  mod <- carehomes$new(p, 0, np, seed = 1L)
  end <- sircovid_date("2020-05-1") / p$dt
  initial <- carehomes_initial(mod$info(), 1, p)
  y <- mod$transform_variables(initial$state)
  i <- c(2, 1, 4, 3)
  y$I_A <- y$I_A[, i, , , drop = FALSE]
  y$T_PCR_pos <- y$T_PCR_pos[, i, , , drop = FALSE]
  y$T_sero_pre_1 <- y$T_sero_pre_1[, i, , , drop = FALSE]
  y$T_sero_pre_2 <- y$T_sero_pre_2[, i, , , drop = FALSE]

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

  z2[["prob_strain"]][, -1] <- z2[["prob_strain"]][2:1, -1, drop = FALSE]
  z2[["cum_sympt_cases_non_variant_over25"]] <-
    z2[["cum_sympt_cases_over25"]] - z2[["cum_sympt_cases_non_variant_over25"]]
  z2$cum_infections_per_strain <-
    z2$cum_infections_per_strain[i, , drop = FALSE]
  ## This one can't easily be computed as it's not quite running
  ## incidence but over a sawtooth; the calculation relative to
  ## cum_infections_per_strain is confirmed elsewhere so here just
  ## move it out the way:
  z2[["sympt_cases_non_variant_over25_inc"]] <-
    z1[["sympt_cases_non_variant_over25_inc"]]
  for (nm in c("T_sero_neg_1", "T_sero_neg_2", "R", "T_PCR_neg")) {
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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  ## Check carehomes_Rt R
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p),
    "Expected R input because there is more than one strain")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[-1, 1, ]),
    "Expected 'R' to have 76 rows")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, -1]),
    "Expected 'R' to have 85 columns")

  ## Check carehomes_Rt prob_strain
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ]),
    "Expected prob_strain input because there is more than one strain")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[-1, 1, ], R = R[, 1, ]),
    "Expected a 2 strains")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ][1, , drop = FALSE],
                 R = R[, 1, ]),
    "Expected 'prob_strain' to have 2 rows")
  expect_error(
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, -1], R = R[, 1, ]),
    "Expected 'prob_strain' to have 85 columns, following 'step'")

  ## Check carehomes_Rt_trajectories R
  expect_error(
    carehomes_Rt_trajectories(steps, S, p),
    "Expected R input because there is more than one strain")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[1, , ], R = R[1, , ]),
    "Expected a 3d array of 'R'")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R[-1, , ]),
    "Expected 'R' to have 76 rows")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R[, -1, ]),
    "Expected 2nd and 3rd")

  ## Check carehomes_Rt_trajectories prob_strain
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, R = R),
    "Expected prob_strain input because there is more than one strain")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[1, , ], R = R),
    "Expected a 3d array of 'prob_strain'")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[-1, , ], R = R),
    "Expected a 3d array of 'prob_strain'")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[, -1, ], R = R),
    "Expected 2nd dim of 'prob_strain' to have length 3, following 'pars'")
  expect_error(
    carehomes_Rt_trajectories(steps, S, p, prob_strain[, , -1], R = R),
    "Expected 3rd dim of 'prob_strain' to have length 85, following 'step'")
})

## Tests for basic object properties, not analytical results from calculations
test_that("wtmean_Rt works as expected", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            cross_immunity = 0)
  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
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

  rt <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ])
  expect_equal(dim(rt$eff_Rt_all), c(85, 2))
  expect_equal(class(rt), c("multi_strain", "Rt"))

  rt_traj <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R)
  expect_equal(dim(rt_traj$eff_Rt_all), c(85, 2, 3))
  expect_equal(class(rt_traj), c("multi_strain", "Rt_trajectories", "Rt"))

  rt_strain_weighted <-
    carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                 weight_Rt = TRUE)
  expect_equal(length(rt_strain_weighted$eff_Rt_all), 85)
  expect_equal(class(rt_strain_weighted), c("single_strain", "Rt"))

  rt_traj_strain_weighted <-
    carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                              weight_Rt = TRUE)
  expect_equal(dim(rt_traj_strain_weighted$eff_Rt_all), c(85, 3))
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
  rt_weight_F <- carehomes_Rt(1, S, p, prob_strain, R = R, weight_Rt = FALSE)
  rt_weight_T <- carehomes_Rt(1, S, p, prob_strain, R = R, weight_Rt = TRUE)
  expect_equal(rt_weight_F$eff_Rt_all[[1]], rt_weight_T$eff_Rt_all)
  expect_approx_equal(rt_weight_F$eff_Rt_general[[1]],
                      rt_weight_T$eff_Rt_general)
  expect_equal(rt_weight_F$Rt_all[[1]], rt_weight_T$Rt_all)
  expect_equal(rt_weight_F$Rt_general[[1]], rt_weight_T$Rt_general)
})


test_that("Can calculate Rt with an empty second variant ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                       weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = FALSE)
  expect_vector_equal(rt_1$eff_Rt_all[, 2], 0)
  expect_vector_equal(rt_1$eff_Rt_general[, 2], 0)
  expect_vector_equal(rt_1$Rt_all[, 2], 0)
  expect_vector_equal(rt_1$Rt_general[, 2], 0)
})




test_that("Can calculate Rt with a second less infectious variant", {
  ## seed with 10 cases on same day as other variant
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0.1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                      weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 5),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)


  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                      weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

  ## Rt should be higher (or equal) for the two variant version
  tol <- 1e-5
  expect_vector_gte(rt_1$Rt_all, rt_1_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_1$Rt_general, rt_1_single_class$Rt_general, tol = tol)
  expect_vector_gte(rt_all$Rt_all, rt_all_single_class$Rt_all, tol = tol)
  expect_vector_gte(rt_all$Rt_general, rt_all_single_class$Rt_general,
                    tol = tol)
})


test_that("Can calculate Rt with a second less letal variant", {
  ## Seed with 10 cases on same day as other variant
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_severity = c(1, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)


  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                       weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = c(1, 0.1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)


  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                       weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_P = c(1, 0.1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)


  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                       weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_C_1 = c(1, 0.1),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)


  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                       weight_Rt = TRUE)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_R))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  R <- y[20:38, , ]

  rt_1_single_class <- carehomes_Rt(steps, S[, 1, ], p, R = R[, 1, ])
  rt_all_single_class <- carehomes_Rt_trajectories(steps, S, p, R = R)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1))

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  info <- mod$info()
  initial <- carehomes_initial(info, 1, p)

  mod$set_state(initial$state, initial$step)
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

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ])
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R)

  ## all values of Rt in 60:70 to be NA and others not to be
  expect_vector_equal(lengths(rt_1[1:3]), 85)
  expect_true(all(is.na(simplify2array(rt_1[4:7])[na_steps, , ])))
  expect_true(!any(is.na(simplify2array(rt_1[4:7])[-na_steps, , ])))

  expect_equal(dim(simplify2array(rt_all[1:3])), c(85, 3, 3))
  expect_true(all(is.na(simplify2array(rt_all[4:7])[na_steps, , , ])))
  expect_true(!any(is.na(simplify2array(rt_all[4:7])[-na_steps, , , ])))


  ## test weighted

  rt_1 <- carehomes_Rt(steps, S[, 1, ], p, prob_strain[, 1, ], R = R[, 1, ],
                       weight_Rt = TRUE)
  rt_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                      weight_Rt = TRUE)

  ## all values of Rt in 60:70 to be NA and others not to be
  expect_vector_equal(lengths(rt_1[1:3]), 85)
  expect_true(all(is.na(simplify2array(rt_1[4:7])[na_steps, ])))
  expect_true(!any(is.na(simplify2array(rt_1[4:7])[-na_steps, ])))

  expect_equal(dim(simplify2array(rt_all[1:3])), c(85, 3, 3))
  expect_true(all(is.na(simplify2array(rt_all[4:7])[na_steps, , ])))
  expect_true(!any(is.na(simplify2array(rt_all[4:7])[-na_steps, , ])))
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

  p <- carehomes_parameters(0, region,
                            waning_rate = 0.1,
                            strain_transmission = c(1, transm_new_variant),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            rel_susceptibility = rep(1, 3),
                            vaccine_progression_rate = numeric(3),
                            vaccine_schedule = vaccine_schedule,
                            vaccine_index_dose2 = 2L,
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 2L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_R <- mod$info()$index$R
  index_prob_strain <- mod$info()$index$prob_strain

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(2)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  R <- y[index_R, , ]
  prob_strain <- y[index_prob_strain, , ]

  for (k in seq_len(np)) { # for each particle

    rt <- carehomes_Rt(steps, S[, k, ], p, prob_strain[, k, ], R = R[, k, ],
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
               "1 or 2")
  expect_error(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_transmission = c(1, 1),
                                    strain_rel_gamma_A = c(2, 5)),
               "must be 1")
  expect_error(carehomes_parameters(sircovid_date("2020-02-07"), "england",
                                    strain_transmission = c(1, 1),
                                    strain_rel_gamma_A = c(1, -1)),
               "non-negative")
})


test_that("Relative gamma = 1 makes no difference", {
  p1 <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  p2 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 1),
                            cross_immunity = 0)
  p3 <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                             strain_transmission = c(1, 1),
                            cross_immunity = 0)
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
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_15_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                         weight_Rt = TRUE)

  ## rate equal to ref
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
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

  rt_1_all <- carehomes_Rt_trajectories(steps, S, p, prob_strain, R = R,
                                        weight_Rt = TRUE)

  ## Rt should be higher (or equal) for the two variant version
  expect_vector_lte(rt_1_all$Rt_all, rt_15_all$Rt_all)
  expect_vector_lte(rt_1_all$Rt_general, rt_15_all$Rt_general)
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
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)
  p$gamma_P_step <- 0

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
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)
  p$gamma_C_1_step <- 0

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
  p$gamma_A_step <- 0
  p$gamma_C_2_step <- 0

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

  ## gammaP is 0 so IC1 is 0 for second strain
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = c(1, 1),
                            strain_rel_gamma_P = c(1, 0),
                            strain_rel_gamma_C_1 = c(1, 1),
                            strain_rel_gamma_C_2 = c(1, 1),
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = c(1, 1),
                            strain_rel_gamma_P = c(1, 1),
                            strain_rel_gamma_C_1 = c(1, 0),
                            strain_rel_gamma_C_2 = c(1, 1),
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, np, seed = 1L)
  mod$set_state(initial$state, initial$step)
  set.seed(1)
  y <- mod$simulate(steps)
  expect_false(all(unlist(y[index_I_C_2_strain_1, , ]) == 0))
  expect_true(all(unlist(y[index_I_C_2_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_1_strain_2, , ]) == 0))


  ## gammaA is 0 & gammaC2 is 0 so R is 0
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_gamma_A = c(1, 0),
                            strain_rel_gamma_P = c(1, 1),
                            strain_rel_gamma_C_1 = c(1, 1),
                            strain_rel_gamma_C_2 = c(1, 0),
                            strain_seed_date =
                              c(sircovid_date("2020-02-07"),
                                sircovid_date("2020-02-08")),
                            strain_seed_rate = c(10, 0),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, np, seed = 1L)
  mod$set_state(initial$state, initial$step)
  set.seed(1)
  y <- mod$simulate(steps)
  expect_false(all(unlist(y[index_R_strain_1, , ]) == 0))
  expect_true(all(unlist(y[index_R_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_2_strain_2, , ]) == 0))
  expect_false(all(unlist(y[index_I_C_2_strain_2, , ]) == 0))
})


test_that("No one is hospitalised, no-one recovers in edge case 2 - multi", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            waning_rate = 1 / 20,
                            strain_seed_rate = c(1, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)
  p$p_C_step[, ] <- 1
  p$p_H_step[, ] <- 1
  p$p_G_D_step[, ] <- 1

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  ## Move initial infectives to 2nd stage sympt
  y0 <- carehomes_initial(info, 1, p)$state
  y0[info$index$I_C_2] <- y0[info$index$I_A]
  y0[info$index$I_A] <- 0

  mod$set_state(y0)
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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(1, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)
  p$p_G_D_step[, ] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$simulate(steps)

  expect_true(all(y[mod$info()$index$G_D, , ] == 0))
})


test_that("G_D strain 2 empty when p_G_D = c(1, 0)", {
  np <- 3L
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_rel_severity = c(1, 0),
                            strain_seed_rate = c(10, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(10, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))
  expect_true(any(y$E[, 1, , , , 2] > 0))
  expect_true(any(y$E[, 2, , , , 2] > 0))
  expect_true(all(y$E[, 3:4, , , , 2] == 0))
})


test_that("Nobody in R2-R4 when strain_transmission = c(1, 0)", {
  np <- 3L
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 0),
                            strain_seed_rate = c(10, 0),
                            beta_value = 1e100,
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_true(all(y$S[, , , 85] == 0))
  expect_true(all(y$R[, 2:4, , , 85] == 0))
  expect_false(all(y$R[, 1, , , 85] == 0))
})


test_that("Can only move to S from R3 and R4 to S", {
  np <- 3L
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            initial_I = 0,
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(10, 0),
                            waning_rate = 1 / 5,
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)
  ## Prevent anyone leaving S
  p$rel_susceptibility[] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)
  info <- mod$info()

  initial <- carehomes_initial(info, np, p)
  y0 <- initial$state
  ## Empty R1 and R2
  y0[info$index$R][1:38] <- 0
  ## Fill R3 and R4
  y0[info$index$R][39:76] <- 1e3

  mod$set_state(y0)

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
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(10, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)
  p$strain_transmission[] <- 1e8
  ## set p_C to 0 so that individuals move to R quickly
  p$p_C_step[, ] <- 0

  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
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
    carehomes_parameters(sircovid_date("2020-02-07"), "england",
                         cross_immunity = 2), "in [0, 1]", fixed = TRUE
  )
  expect_error(
    carehomes_parameters(sircovid_date("2020-02-07"), "england",
                         cross_immunity = -2), "in [0, 1]", fixed = TRUE
  )
  expect_error(
    carehomes_parameters(sircovid_date("2020-02-07"), "england",
                         cross_immunity = c(1, 1)), "Invalid length"
  )
})

test_that("complete cross_immunity means no Strain 3/4 infections", {
  np <- 1L

  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(10, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 1)
  mod <- carehomes$new(p, 0, np)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_true(all(y$cum_infections_per_strain[3:4, ] == 0))
})

test_that("some cross-immunity means less Strain 3 or 4 infections than none
           and > 0", {
  np <- 1L

  ## no cross-immnunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(10, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0)
  mod <- carehomes$new(p, 0, np)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))
  infect_no_cross <- y$cum_infections_per_strain[3:4, 85]

  ## some cross-immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(10, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = 0.5)
  mod <- carehomes$new(p, 0, np)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))
  infect_some_cross <- y$cum_infections_per_strain[3:4, 85]

  expect_vector_lte(infect_some_cross, infect_no_cross)
  expect_true(all(infect_some_cross > 0))
  expect_true(all(infect_no_cross > 0))
})


test_that("cross-immunity can be separated by strain", {
  seed <- 1
  np <- 1L
  set.seed(seed)

  ## complete immunity from Strain 1 means Strain 3 empty (1 -> 2)
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
    strain_transmission = c(1, 1),
    strain_seed_rate = c(10, 0),
    strain_seed_date =
      sircovid_date(c("2020-02-07", "2020-02-08")),
    cross_immunity = c(1, 0)
  )
  mod <- carehomes$new(p, 0, np, seed = seed)
  initial <- carehomes_initial(mod$info(), np, p)
  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)
  mod$set_state(initial$state, initial$step)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_equal(y$cum_infections_per_strain[3, 85], 0)
  expect_gt(y$cum_infections_per_strain[4, 85], 0)

  ## complete immunity from Strain 2 means Strain 4 empty (2 -> 1)
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            strain_seed_rate = c(100, 0),
                            strain_seed_date =
                              sircovid_date(c("2020-02-07", "2020-02-08")),
                            cross_immunity = c(0, 1))
  mod <- carehomes$new(p, 0, np, seed = seed)
  initial <- carehomes_initial(mod$info(), np, p)
  mod$set_state(initial$state, initial$step)
  set.seed(1)
  y <- mod$transform_variables(
    drop(mod$simulate(steps)))

  expect_equal(y$cum_infections_per_strain[4, 85], 0)
  expect_gt(y$cum_infections_per_strain[3, 85], 0)
})


test_that("Can calculate ifr_t with an empty second variant ", {
  ## Run model with 2 variants, but both have same transmissibility
  ## no seeding for second variant so noone infected with that one
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            strain_transmission = c(1, 1),
                            rel_susceptibility = c(1, 1, 1),
                            cross_immunity = 0)

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  I <- y[index_I, , ]

  ifr_t_1 <- carehomes_ifr_t(steps, S[, 1, ], I[, 1, ], p)
  ifr_t_all <- carehomes_ifr_t_trajectories(steps, S, I, p)

  ## Run model with one strain only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")

  np <- 3L
  mod <- carehomes$new(p, 0, np, seed = 1L)

  initial <- carehomes_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted

  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  mod$set_index(c(index_S, index_I))
  y <- mod$simulate(steps)

  S <- y[1:19, , ]
  I <- y[20:38, , ]

  ifr_t_1_single_class <- carehomes_ifr_t(steps, S[, 1, ], I[, 1, ], p)
  ifr_t_all_single_class <- carehomes_ifr_t_trajectories(steps, S, I, p)

  expect_equal(ifr_t_1$step, ifr_t_1_single_class$step)
  expect_equal(ifr_t_1$date, ifr_t_1_single_class$date)
  expect_equal(ifr_t_1$beta, ifr_t_1_single_class$beta)
  expect_equal(ifr_t_1$eff_ifr_t_all, ifr_t_1_single_class$eff_ifr_t_all)
  expect_equal(ifr_t_1$eff_Rt_general,
               ifr_t_1_single_class$eff_Rt_general)
  expect_equal(ifr_t_1$ifr_t_all, ifr_t_1_single_class$ifr_t_all)
  expect_equal(ifr_t_1$Rt_general, ifr_t_1_single_class$Rt_general)

  expect_equal(ifr_t_all$step, ifr_t_all_single_class$step)
  expect_equal(ifr_t_all$date, ifr_t_all_single_class$date)
  expect_equal(ifr_t_all$beta, ifr_t_all_single_class$beta)
  expect_equal(ifr_t_all$eff_ifr_t_all,
               ifr_t_all_single_class$eff_ifr_t_all)
  expect_equal(ifr_t_all$eff_Rt_general,
               ifr_t_all_single_class$eff_Rt_general)
  expect_equal(ifr_t_all$ifr_t_all, ifr_t_all_single_class$ifr_t_all)
  expect_equal(ifr_t_all$Rt_general,
               ifr_t_all_single_class$Rt_general)

})

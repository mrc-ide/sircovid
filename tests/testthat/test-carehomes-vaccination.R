context("carehomes (vaccination)")

test_that("there are no infections if everyone is vaccinated with a vaccine
          preventing 100% of acquisition", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", rel_susceptibility = c(1, 0),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  mod$set_state(state)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## Reshape to show the full shape of s
  expect_equal(length(s), prod(info$dim$S) * 101)
  s <- array(s, c(info$dim$S, 101))

  ## Noone moves into unvaccinated
  ## except in the group where infections because of waning immunity
  expect_true(all(s[-4, 1, ] == 0))

  ## Noone changes compartment within the vaccinated individuals
  expect_true(all(s[, 2, ] == s[, 2, 1]))
})

test_that("Everyone moves to vaccinated and stays there if everyone quickly
          gets vaccinated with a vaccine preventing 100% of acquisition
          and no waning immunity", {
            p <- carehomes_parameters(0, "england",
                                      beta_value = c(0, 0, 1),
                                      beta_date = c(0, 4, 5),
                                      rel_susceptibility = c(1, 0),
                                      vaccine_progression_rate = c(Inf, 0))
            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()
            mod$set_state(carehomes_initial(info, 1, p)$state)
            mod$set_index(integer(0))
            y <- mod$transform_variables(drop(
              dust::dust_iterate(mod, seq(0, 400, by = 4))))
            expect_true(all(y$S[, 1, 1] == y$S[, 2, 2] + rowSums(y$E[, , , 2])))
            expect_true(all(y$S[, , 101] == y$S[, , 2]))
})


test_that("Everyone moves to waning immunity stage and stays there if everyone
          quickly gets vaccinated and loses immunity", {
            p <- carehomes_parameters(0, "england",
                                      beta_value = 0,
                                      rel_susceptibility = c(1, 0, 0),
                                      vaccine_progression_rate = c(Inf, Inf, 0))
            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()
            mod$set_state(carehomes_initial(info, 1, p)$state)
            mod$set_index(integer(0))
            y <- mod$transform_variables(drop(
              dust::dust_iterate(mod, seq(0, 400, by = 4))))
            expect_true(all(y$S[, 1, 1] == y$S[, 3, 2]))
            expect_true(all(y$S[, , 101] == y$S[, , 2]))
})

test_that("Everyone moves to last of 5 waning immunity stage and stays there if
          everyone quickly gets vaccinated and loses immunity", {
            p <- carehomes_parameters(0, "england",
                                      beta_value = 0,
                                      rel_susceptibility = c(1, rep(0, 10), 0),
                                      vaccine_progression_rate =
                                        c(Inf, rep(Inf, 10), 0))
            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()
            mod$set_state(carehomes_initial(info, 1, p)$state)
            mod$set_index(integer(0))
            y <- mod$transform_variables(drop(
              dust::dust_iterate(mod, seq(0, 400, by = 12))))
            expect_true(all(y$S[, 1, 1] == y$S[, 12, 2]))
            expect_true(all(y$S[, , 34] == y$S[, , 2]))
})

test_that("Everyone moves to the R compartment corresponding to their
          vaccination class if there is a high beta and no vaccine
          progression", {
            p <- carehomes_parameters(0, "england",
                                      beta_value = 1e9,
                                      rel_susceptibility = c(1, 1, 1),
                                      vaccine_progression_rate =
                                        c(0, 0, 0))

            # increase progression rates
            p[grep("gamma", names(p))] <- 1e9
            # make p_death zero
            param_death_idx <-
              intersect(grep("p_death", names(p)), grep("_step", names(p)))
            p[param_death_idx] <- 0

            mod <- carehomes$new(p, 0, 1)
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
            mod$set_index(integer(0))
            y <- dust::dust_iterate(mod, seq(0, 400, by = 4), index)

            ## Reshape to show the full shape of s
            expect_equal(length(y),
                         prod(info$dim$S) * 101 + prod(info$dim$R) * 101)
            s <- array(y[seq_len(prod(info$dim$S)), , ], c(info$dim$S, 101))
            r <- array(y[prod(info$dim$S) + seq_len(prod(info$dim$R)), , ],
                       c(info$dim$R, 101))

            ## all have moved from S to R in relevant vaccination class
            ## ignoring age group 4 where infections are seeded
            expect_true(all(s[-4, , 1] == r[-4, , 101]))
          })

test_that("Everyone moves back to unvaccinated from vacinated if large waning
          of immunity and no vaccination", {
            p <- carehomes_parameters(0, "england",
                                      beta_value = 0,
                                      rel_susceptibility = c(1, 0),
                                      vaccine_progression_rate = c(0, Inf))

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_S <- array(info$index$S, info$dim$S)
            state[index_S[, 2]] <- state[index_S[, 1]]
            state[index_S[, 1]] <- 0

            mod$set_state(state)
            mod$set_index(integer(0))
            s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

            ## Reshape to show the full shape of s
            expect_equal(length(s), prod(info$dim$S) * 101)
            s <- array(s, c(info$dim$S, 101))

            ## noone left in vaccinated at time step 2
            expect_true(all(s[, 2, 2] == 0))

            ## everybody back in unvaccinated at time step 2
            expect_true(all(s[, 1, 2] == s[, 2, 1]))
})

test_that("there are no vaccinations when vaccination rate is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  waning_rate <- rep(1 / 20, 19)
  waning_rate[4] <- 0 # no waning in group with seeded infections
  # otherwise S can go up as these infected individuals loose immunity

  p <- carehomes_parameters(0, "england",
                            waning_rate = waning_rate)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## No vaccinated susceptibles:
  expect_equal(s[-seq_len(carehomes_n_groups()), , ],
               array(0, c(nrow(s) - carehomes_n_groups(), 101)))
})

test_that("Can calculate Rt with an (empty) vaccination class", {

  ## run model with unvaccinated & vaccinated, but both have same susceptibility
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1))

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

  ## run model with unvaccinated class only
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, 1))

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

test_that("Effective Rt reduced by rel_susceptbility if all vaccinated", {

  reduced_susceptibility <- 0.2 # can put anything here

  ## run model with unvaccinated & vaccinated (with susceptibility halved)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england",
                            rel_susceptibility = c(1, reduced_susceptibility),
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
  y <- dust::dust_iterate(mod, steps, index)

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


test_that("N_tot, N_tot2 and N_tot3 stay constant with vaccination", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  set.seed(1)
  p <- carehomes_parameters(0, "uk", waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5, 0.1),
                            vaccine_progression_rate = c(1, 0.5, 0.01))

  p$p_asympt <- rep(1, 19)

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

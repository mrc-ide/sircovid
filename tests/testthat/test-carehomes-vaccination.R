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
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  
  state <- carehomes_initial(info, 1, p)$state
  
  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0
  
  mod$set_state(state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(
    dust::dust_iterate(mod, seq(0, 400, by = 4))))
  
  ## Noone moves into I_ILI or I_mild ever
  ## other than in the 4th age group where some infections are seeded
  ## in the unvaccinated group and because of waning immunity they may 
  ## eventually end up in I_ILI or I_mild upon reinfection
  expect_true(all(y$I_ILI[-4, , , ] == 0))
  expect_true(all(y$I_mild[-4, , , ] == 0))
  
})


test_that("Noone hospitalised with perfect vaccine wrt rel_p_hosp_if_sympt_if_sympt", {
  ## i.e. if everyone is vaccinated with a vaccine preventing
  ## 100% of hospitalisations
  
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", 
                            rel_susceptibility = c(1, 1),
                            rel_p_sympt = c(1, 1),
                            rel_p_hosp_if_sympt = c(1, 0),
                            waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  
  state <- carehomes_initial(info, 1, p)$state
  
  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0
  
  mod$set_state(state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(
    dust::dust_iterate(mod, seq(0, 400, by = 4))))
  
  ## Noone moves into hospitalised compartments ever
  ## other than in the 4th age group where some infections are seeded
  ## in the unvaccinated group and because of waning immunity they may 
  ## eventually end up in hospital upon reinfection
  expect_true(all(y$D_hosp[-4, ] == 0))
  expect_true(all(y$I_hosp_R_unconf[-4, , , ] == 0))
  expect_true(all(y$I_hosp_R_conf[-4, , , ] == 0))
  expect_true(all(y$I_hosp_D_unconf[-4, , , ] == 0))
  expect_true(all(y$I_hosp_D_conf[-4, , , ] == 0))
})


test_that("Vaccination of susceptibles works", {
  ## Tests that:
  ## Every susceptible moves to vaccinated and stays there if
  ## everyone quickly gets vaccinated with a vaccine preventing 100% of
  ## acquisition and no waning immunity
            p <- carehomes_parameters(0, "england",
                                      beta_value = c(0, 0, 1),
                                      beta_date = c(0, 4, 5),
                                      rel_susceptibility = c(1, 0),
                                      rel_p_sympt = c(1, 1),
                                      rel_p_hosp_if_sympt = c(1, 1),
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


test_that("Vaccination of exposed individuals works", {
  ## Tests that:
  ## Every exposed moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if disease progression is stopped after E
            p <- carehomes_parameters(0, "england",
                                      rel_susceptibility = c(1, 0),
                                      rel_p_sympt = c(1, 1),
                                      rel_p_hosp_if_sympt = c(1, 1),
                                      vaccine_progression_rate = c(Inf, 0))

            # stop disease progression after E
            p$gamma_E <- 0

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_E <- array(info$index$E, info$dim$E)
            index_S <- array(info$index$S, info$dim$S)
            state[index_E[, 1, ]] <- state[index_S]
            state[index_E[, 2, ]] <- state[index_S]
            state[index_S] <- 0

            mod$set_state(state)
            mod$set_index(integer(0))
            e <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$E)

            ## Reshape to show the full shape of e
            expect_equal(length(e), prod(info$dim$E) * 101)
            e <- array(e, c(info$dim$E, 101))

            ## every E moves from unvaccinated to vaccinated between
            ## time steps 1 and 2
            E_compartment_idx <- 1
            unvacc_idx <- 1
            vacc_idx <- 2
            expect_true(all(e[, E_compartment_idx, unvacc_idx, 1] ==
                              e[, E_compartment_idx, vacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(e[, E_compartment_idx, vacc_idx, 2] ==
                              e[, E_compartment_idx, vacc_idx, 101]))
            ## same for second E compartment
            E_compartment_idx <- 2
            expect_true(all(e[, E_compartment_idx, unvacc_idx, 1] ==
                              e[, E_compartment_idx, vacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(e[, E_compartment_idx, vacc_idx, 2] ==
                              e[, E_compartment_idx, vacc_idx, 101]))
})


test_that("Vaccination of asymptomatic infectious individuals works", {
  ## Tests that:
  ## Every I_asympt moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if disease progression is stopped after I_asympt
            p <- carehomes_parameters(0, "england",
                                      rel_susceptibility = c(1, 0),
                                      rel_p_sympt = c(1, 1),
                                      rel_p_hosp_if_sympt = c(1, 1),
                                      vaccine_progression_rate = c(Inf, 0))

            # stop disease progression after I_asympt
            p$gamma_asympt <- 0

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_I_asympt <- array(info$index$I_asympt, info$dim$I_asympt)
            index_S <- array(info$index$S, info$dim$S)
            state[index_I_asympt[, 1, ]] <- state[index_S]
            state[index_S] <- 0

            mod$set_state(state)
            mod$set_index(integer(0))
            i_asympt <-
              dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$I_asympt)

            ## Reshape to show the full shape of i_asympt
            expect_equal(length(i_asympt), prod(info$dim$I_asympt) * 101)
            i_asympt <- array(i_asympt, c(info$dim$I_asympt, 101))

            ## every I_asympt moves from unvaccinated to vaccinated between
            ## time steps 1 and 2
            I_asympt_compartment_idx <- 1
            unvacc_idx <- 1
            vacc_idx <- 2
            expect_true(all(
              i_asympt[, I_asympt_compartment_idx, unvacc_idx, 1] ==
                i_asympt[, I_asympt_compartment_idx, vacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(
              i_asympt[, I_asympt_compartment_idx, vacc_idx, 2] ==
                i_asympt[, I_asympt_compartment_idx, vacc_idx, 101]))
})


test_that("Vaccination of recovered individuals works", {
  ## Test that:
  ## Every R moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a vaccine with no waning immunity
  ## and if no natural waning of immunity
            p <- carehomes_parameters(0, "england",
                                      rel_susceptibility = c(1, 0),
                                      rel_p_sympt = c(1, 1),
                                      rel_p_hosp_if_sympt = c(1, 1),
                                      vaccine_progression_rate = c(Inf, 0))

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_I_asympt <- array(info$index$I_asympt, info$dim$I_asympt)
            index_R <- array(info$index$R, info$dim$R)
            index_R_neg <- array(info$index$R_neg, info$dim$R_neg)
            index_PCR_neg <- array(info$index$PCR_neg, info$dim$PCR_neg)
            index_S <- array(info$index$S, info$dim$S)
            state[index_R] <- state[index_S]
            state[index_R_neg] <- state[index_S]
            state[index_PCR_neg] <- state[index_S]
            state[index_S] <- 0
            state[index_I_asympt] <- 0 # remove seeded infections

            mod$set_state(state)
            mod$set_index(integer(0))
            r <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$R)

            ## Reshape to show the full shape of r
            expect_equal(length(r), prod(info$dim$R) * 101)
            r <- array(r, c(info$dim$R, 101))

            ## every R moves from unvaccinated to vaccinated between
            ## time steps 1 and 2
            unvacc_idx <- 1
            vacc_idx <- 2
            expect_true(all(r[, unvacc_idx, 1] == r[, vacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(r[, vacc_idx, 2] == r[, vacc_idx, 101]))
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

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_E <- array(info$index$E, info$dim$E)
            index_S <- array(info$index$S, info$dim$S)
            state[index_E[, 1, 2]] <- state[index_S[, 1]]
            state[index_E[, 1, 1]] <- 0
            state[index_E[, 2, 2]] <- state[index_S[, 1]]
            state[index_E[, 2, 1]] <- 0
            state[index_S] <- 0

            mod$set_state(state)
            mod$set_index(integer(0))
            e <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$E)

            ## Reshape to show the full shape of e
            expect_equal(length(e), prod(info$dim$E) * 101)
            e <- array(e, c(info$dim$E, 101))

            ## every E moves from vaccinated to unvaccinated between
            ## time steps 1 and 2
            E_compartment_idx <- 1
            unvacc_idx <- 1
            vacc_idx <- 2
            expect_true(all(e[, E_compartment_idx, vacc_idx, 1] ==
                              e[, E_compartment_idx, unvacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(e[, E_compartment_idx, unvacc_idx, 2] ==
                              e[, E_compartment_idx, unvacc_idx, 101]))
            ## same for second E compartment
            E_compartment_idx <- 2
            expect_true(all(e[, E_compartment_idx, vacc_idx, 1] ==
                              e[, E_compartment_idx, unvacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(e[, E_compartment_idx, unvacc_idx, 2] ==
                              e[, E_compartment_idx, unvacc_idx, 101]))
})


test_that("Returning to unvaccinated stage works for I_asympt individuals", {
  ## Tests that:
  ## Every I_asympt moves back from vaccinated to unvaccinated and stays
  ## there if vaccine has fast waning immunity,
  ## beta is zero, and if disease progression
  ## is stopped after I_asympt
            p <- carehomes_parameters(0, "england",
                                      beta_value = 0,
                                      rel_susceptibility = c(1, 0),
                                      rel_p_sympt = c(1, 1),
                                      rel_p_hosp_if_sympt = c(1, 1),
                                      vaccine_progression_rate = c(0, Inf))

            # stop disease progression after I_asympt
            p$gamma_asympt <- 0

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_I_asympt <- array(info$index$I_asympt, info$dim$I_asympt)
            index_S <- array(info$index$S, info$dim$S)
            state[index_I_asympt[, 1, 2]] <- state[index_S[, 1]]
            state[index_I_asympt[, 1, 1]] <- 0
            state[index_S] <- 0

            mod$set_state(state)
            mod$set_index(integer(0))
            i_asympt <- dust::dust_iterate(mod, seq(0, 400, by = 4),
                                           info$index$I_asympt)

            ## Reshape to show the full shape of i_asympt
            expect_equal(length(i_asympt), prod(info$dim$I_asympt) * 101)
            i_asympt <- array(i_asympt, c(info$dim$I_asympt, 101))

            ## every I_asympt moves from vaccinated to unvaccinated between
            ## time steps 1 and 2
            I_asympt_compartment_idx <- 1
            unvacc_idx <- 1
            vacc_idx <- 2
            expect_true(all(
              i_asympt[, I_asympt_compartment_idx, vacc_idx, 1] ==
                i_asympt[, I_asympt_compartment_idx, unvacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(
              i_asympt[, I_asympt_compartment_idx, unvacc_idx, 2] ==
                i_asympt[, I_asympt_compartment_idx, unvacc_idx, 101]))
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

            mod <- carehomes$new(p, 0, 1)
            info <- mod$info()

            state <- carehomes_initial(info, 1, p)$state

            index_I_asympt <- array(info$index$I_asympt, info$dim$I_asympt)
            index_R <- array(info$index$R, info$dim$R)
            index_R_neg <- array(info$index$R_neg, info$dim$R_neg)
            index_PCR_neg <- array(info$index$PCR_neg, info$dim$PCR_neg)
            index_S <- array(info$index$S, info$dim$S)
            state[index_R[, 2]] <- state[index_S[, 1]]
            state[index_R[, 1]] <- 0
            state[index_R_neg[, 2]] <- state[index_S[, 1]]
            state[index_R_neg[, 1]] <- 0
            state[index_PCR_neg[, 2]] <- state[index_S[, 1]]
            state[index_PCR_neg[, 1]] <- 0
            state[index_S] <- 0
            state[index_I_asympt] <- 0 # remove seeded infections

            mod$set_state(state)
            mod$set_index(integer(0))
            r <- dust::dust_iterate(mod, seq(0, 400, by = 4),
                                           info$index$R)

            ## Reshape to show the full shape of r
            expect_equal(length(r), prod(info$dim$R) * 101)
            r <- array(r, c(info$dim$R, 101))

            ## every R moves from vaccinated to unvaccinated between
            ## time steps 1 and 2
            unvacc_idx <- 1
            vacc_idx <- 2
            expect_true(all(r[, vacc_idx, 1] == r[, unvacc_idx, 2]))
            ## then they don't move anymore
            expect_true(all(r[, unvacc_idx, 2] == r[, unvacc_idx, 101]))
})


test_that("Vaccine progression through 3 classes works for susceptibles", {
  ## Tests that:
  ## Every susceptible moves to waning immunity stage and stays there if
  ## everyone quickly gets vaccinated and loses immunity
            p <- carehomes_parameters(0, "england",
                                      beta_value = 0,
                                      rel_susceptibility = c(1, 0, 0),
                                      rel_p_sympt = c(1, 1, 1),
                                      rel_p_hosp_if_sympt = c(1, 1, 1),
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

test_that("Vaccine progression through 12 classes works for susceptibles", {
  ## Tests that:
  ## Every susceptible moves to last of 12 waning immunity stage and stays
  ## there if everyone quickly gets vaccinated and loses immunity
            p <- carehomes_parameters(0, "england",
                                      beta_value = 0,
                                      rel_susceptibility = c(1, rep(0, 10), 0),
                                      rel_p_sympt = rep(1, 12),
                                      rel_p_hosp_if_sympt = rep(1, 12),
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

test_that("Clinical progression within a vaccination class works", {
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

test_that("there are no vaccinated susceptibles when vaccination rate is 0", {
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
  y <- dust::dust_iterate(mod, steps, index)

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
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(1, 0.5, 0.01))

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


test_that(
  "N_tot, N_tot2 and N_tot3 stay constant with high rates of vaccination", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  set.seed(1)
  ## TODO: set up a more specific set of tests to test the combined moves
  ## whereby in a single times step an individual progresses to next clinical
  ## stage and progresses to the next vaccination stage
  p <- carehomes_parameters(0, "uk", waning_rate = 1 / 20,
                            rel_susceptibility = c(1, 0.5, 0.1),
                            rel_p_sympt = c(1, 1, 1),
                            rel_p_hosp_if_sympt = c(1, 1, 1),
                            vaccine_progression_rate = c(500, 100, 50))

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

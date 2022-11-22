context("lancelot (check)")

test_that("N_tots stay constant without waning immunity", {
  ## waning_rate default is 0
  p <- lancelot_parameters(0, "uk")
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)
  mod$update_state(state = lancelot_initial(info, 1, p))
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

test_that("N_tot stays constant with waning immuity, while sero and PCR N_tots
          are non-decreasing", {
            p <- lancelot_parameters(0, "uk", waning_rate = 1 / 20)
            mod <- lancelot$new(p, 0, 1)
            info <- mod$info()
            y0 <- lancelot_initial(info, 1, p)
            mod$update_state(state = lancelot_initial(info, 1, p))
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

test_that("N_tot stays constant when p_R < 1, while
          sero and PCR N_tots are non-decreasing", {
            p <- lancelot_parameters(0, "uk")
            p$p_R_step[, ] <- 0.5
            mod <- lancelot$new(p, 0, 1)
            info <- mod$info()
            y0 <- lancelot_initial(info, 1, p)
            mod$update_state(state = lancelot_initial(info, 1, p))
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

test_that("there are no infections when beta is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  waning_rate <- rep(1 / 20, 19)
  waning_rate[4] <- 0 # no waning in group with seeded infections
  # otherwise S can go up as these infected individuals loose immunity

  p <- lancelot_parameters(0, "england", beta_value = 0,
                           waning_rate = waning_rate)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Susceptible population is never drawn down after initial seeding:
  expect_equal(s[, , -1, drop = FALSE], array(s[, , 2], c(19, 1, 100)))
})


test_that("everyone is infected when beta is large", {
  p <- lancelot_parameters(0, "england", beta_value = 1e9)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  state <- lancelot_initial(info, 1, p)

  ## seed directly into I_A
  index_I_A <- array(info$index$I_A, info$dim$I_A)
  state[index_I_A] <- 5

  mod$update_state(state = state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  expect_vector_equal(s[, , -1], 0)
})


test_that("noone stays in R if waning rate is very
          large", {
  # with a large waning rate and beta = 0,
  # people can move from R to S but not outside of S
  # therefore R should quickly get empty
  p <- lancelot_parameters(0, "england",
                           beta_value = 0, # to forbid movement out of S
                           waning_rate = Inf) # to force movement out of R
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  state <- lancelot_initial(info, 1, p)

  # move everyone to R
  index_S <- array(info$index$S, info$dim$S)
  index_R <- array(info$index$R, info$dim$R)
  state[index_R] <- rowSums(array(state[index_S], info$dim$S))
  state[index_S] <- 0

  mod$update_state(state = state)
  y <- mod$transform_variables(drop(
    mod$simulate(seq(0, 400, by = 4))))

  # other than in the 4th age group (where infections are seeded)
  # after the first day (4 times steps), R is empty
  expect_true(all(y$R[-4, , -1, ] == 0))

})

test_that("R is non-decreasing and S is non-increasing if waning rate is 0", {
  p <- lancelot_parameters(0, "england",
                           waning_rate = 0)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(drop(
    mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(diff(t(drop(y$R))) >= 0))
  expect_true(all(diff(t(y$S[, 1, ])) <= 0))

})


test_that("No one is infected if there is no seeding", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20,
                           initial_seed_size = 0)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y <- lancelot_initial(info, 1, p)
  y[info$index$I_A] <- 0
  y[info$index$T_sero_pre_1] <- 0
  y[info$index$T_sero_pre_2] <- 0
  y[info$index$T_PCR_pos] <- 0

  mod$update_state(state = y)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(nrow(s), 1, 101)))
})


test_that("No one is hospitalised, no-one dies if p_C is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_C_step[, ] <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$E > 0L))
  expect_true(all(y$I_P == 0))
  expect_true(all(y$I_C_1 == 0))
  expect_true(all(y$I_C_2 == 0))
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
  expect_true(all(y$D_hosp == 0))
  expect_true(all(y$G_D == 0))
  expect_true(all(y$D_non_hosp == 0))
})


test_that("No one is hospitalised, no-one dies if p_H is 0", {
  set.seed(1)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_H_step[, ] <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$E > 0L))
  expect_true(any(y$I_P > 0))
  expect_true(any(y$I_C_1 > 0))
  expect_true(any(y$I_C_2 > 0))
  expect_true(all(y$H_R_unconf == 0))
  expect_true(all(y$H_R_conf == 0))
  expect_true(all(y$H_D_unconf == 0))
  expect_true(all(y$H_D_conf == 0))
  expect_true(all(y$I_ICU_R_unconf == 0))
  expect_true(all(y$I_ICU_R_conf == 0))
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
  expect_true(all(y$D_hosp == 0))
  expect_true(all(y$G_D == 0))
  expect_true(all(y$D_non_hosp == 0))
})


test_that("No one is hospitalised, no-one recovers in edge case", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_H_step[, ] <- 1
  p$p_G_D_step[, ] <- 1
  p$p_C_step[, ] <- 1

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  ## Move initial infectives to 2nd stage sympt
  y0 <- lancelot_initial(info, 1, p)
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


test_that("No one is hospitalised, no-one recovers in edge case 2", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_H_step[, ] <- 1
  p$p_G_D_step[, ] <- 1
  p$p_C_step[, ] <- 1

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  ## Move initial infectives to 2nd stage sympt
  y0 <- lancelot_initial(info, 1, p)
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

test_that("No-one recovers if p_R = 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_R_step[, ] <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  ## Add initial individuals into compartments that feed into R
  y0 <- lancelot_initial(info, 1, p)
  y0[info$index$I_A] <- 50
  y0[info$index$I_C_2] <- 50
  y0[info$index$H_R_unconf] <- 50
  y0[info$index$H_R_conf] <- 50
  y0[info$index$W_R_unconf] <- 50
  y0[info$index$W_R_conf] <- 50

  mod$update_state(state = y0)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_A > 0))
  expect_true(any(y$I_C_2 > 0))
  expect_true(any(y$H_R_unconf > 0))
  expect_true(any(y$H_R_conf > 0))
  expect_true(any(y$W_R_unconf > 0))
  expect_true(any(y$W_R_conf > 0))
  expect_true(all(y$R == 0))
})

test_that("Everyone recovers in edge case", {
  ## This test is primarily to test the behaviour for p_R = 1
  p <- lancelot_parameters(0, "england", beta_value = 1e9)
  p$p_R_step[, ] <- 1
  p$p_G_D_step[, ] <- 0
  p$p_ICU_D_step[, ] <- 0
  p$p_H_D_step[, ] <- 0
  p$p_W_D_step[, ] <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()

  ## Add initial individuals into compartments that feed into R
  y0 <- lancelot_initial(info, 1, p)
  y0[info$index$I_A] <- 50
  y0[info$index$I_C_2] <- 50
  y0[info$index$H_R_unconf] <- 50
  y0[info$index$H_R_conf] <- 50
  y0[info$index$W_R_unconf] <- 50
  y0[info$index$W_R_conf] <- 50

  mod$update_state(state = y0)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_A > 0))
  expect_true(any(y$I_C_2 > 0))
  expect_true(any(y$H_R_unconf > 0))
  expect_true(any(y$H_R_conf > 0))
  expect_true(any(y$W_R_unconf > 0))
  expect_true(any(y$W_R_conf > 0))
  expect_true(all(y$H_D_unconf == 0))
  expect_true(all(y$H_D_conf == 0))
  expect_true(all(y$W_D_unconf == 0))
  expect_true(all(y$W_D_conf == 0))
  expect_true(all(y$W_D_unconf == 0))
  expect_true(all(y$ICU_D_conf == 0))
  expect_true(all(y$ICU_D_unconf == 0))
  expect_true(all(y$G_D == 0))
  expect_true(any(y$R > 0))

  expect_true(all(apply(y$S, c(1, 2), diff) <= 0))
})

test_that("No one dies in the community if p_G_D is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_G_D_step[, ] <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_P + y$I_A > 0))
  expect_true(all(y$G_D == 0))
  expect_true(all(y$D_non_hosp == 0))
})


test_that("forcing hospital route results in correct path", {
  ## We're going to try a number of very similar tests here, so a
  ## helper function will help run a model with given probabilities
  ## and verify that some compartments have cases (nonzero) and others
  ## area all zeros.
  helper <- function(prob_ICU, prob_ICU_D, prob_H_D,
                     prob_W_D, expect_cases, expect_zero) {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity
    p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
    p$p_ICU_step[, ] <- prob_ICU %||% p$p_ICU_step
    p$p_ICU_D_step[, ] <- prob_ICU_D %||% p$p_ICU_D_step
    p$p_H_D_step[, ] <- prob_H_D %||% p$p_H_D_step
    p$p_W_D_step[, ] <- prob_W_D %||% p$p_W_D_step

    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    mod$update_state(state = lancelot_initial(info, 1, p))
    y <- mod$transform_variables(
      drop(mod$simulate(seq(0, 400, by = 4))))

    ## Save some work by using the total of confirmed and unconfirmed
    y$H_R <- y$H_R_unconf + y$H_R_conf
    y$H_D <- y$H_D_unconf + y$H_D_conf
    y$ICU_pre <- y$ICU_pre_unconf + y$ICU_pre_conf
    y$ICU_W_R <- y$ICU_W_R_unconf + y$ICU_W_R_conf
    y$ICU_W_D <- y$ICU_W_D_unconf + y$ICU_W_D_conf
    y$ICU_D <- y$ICU_D_unconf + y$ICU_D_conf
    y$W_R <- y$W_R_unconf + y$W_R_conf
    y$W_D <- y$W_D_unconf + y$W_D_conf

    for (i in expect_cases) {
      expect_true(any(y[[i]] > 0), label = sprintf("Expected cases in %s", i))
    }
    for (i in expect_zero) {
      expect_true(all(y[[i]] == 0), label = sprintf("Expected zeros in %s", i))
    }
  }

  ## p_ICU = 0, p_H_D = 0 no-one goes into ICU, no deaths
  helper(0, NULL, 0, NULL, "H_R",
         c("H_D", "ICU_W_R", "ICU_W_D", "ICU_D", "ICU_pre",
           "W_R", "W_D", "D_hosp"))

  ## p_death_hosp = 1, p_ICU = 0 no-one goes into ICU, no
  ## recovery in hospital
  helper(0, NULL, 1, NULL, "H_D",
         c("H_R", "ICU_W_R", "ICU_W_D", "ICU_D", "ICU_pre",
           "W_R", "W_D"))

  ## p_ICU_D = 1, p_ICU = 1 no-one goes in hosp_D / hosp_R,
  ## no recovery from ICU
  helper(1, 1, NULL, NULL, "ICU_D",
         c("H_R", "H_D", "ICU_W_R", "ICU_W_D", "W_R",
           "W_D"))

  ## p_ICU_D = 0, p_ICU = 1, p_W_D = 0 no-one goes in
  ## hosp_D / hosp_R, no deaths
  helper(1, 0, NULL, 0, c("ICU_W_R", "W_R"),
         c("H_R", "H_D", "ICU_W_D", "ICU_D", "W_D",
           "D_hosp"))

  ## p_ICU_D = 0, p_ICU = 1, p_W_D = 1 no-one goes in
  ## hosp_D / hosp_R, no-one recovers in stepdown
  helper(1, 0, NULL, 1, c("ICU_W_D", "W_D"),
         c("H_R", "H_D", "ICU_W_R", "ICU_D", "W_R"))
})


test_that("No one seroconverts if p_sero_pos is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_sero_pos_1[] <- 0
  p$p_sero_pos_2[] <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$T_sero_neg_1 > 0))
  expect_true(all(y$T_sero_pos_1 == 0))
  expect_true(any(y$T_sero_neg_2 > 0))
  expect_true(all(y$T_sero_pos_2 == 0))
})


test_that("No one does not seroconvert and no one seroreverts
          if p_sero_pos is 1 and gamma_sero_pos is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)

  p$p_sero_pos_1[] <- 1
  p$p_sero_pos_2[] <- 1
  ## set gamma_sero_pos_1 to 0 so no-one seroreverts
  p$gamma_sero_pos_1 <- 0
  p$gamma_sero_pos_2 <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$T_sero_neg_1 == 0))
  expect_true(any(y$T_sero_pos_1 > 0))
  expect_true(all(y$T_sero_neg_2 == 0))
  expect_true(any(y$T_sero_pos_2 > 0))
})


test_that("setting a gamma to Inf results immediate progression", {
  helper <- function(gamma_name, progression_name, compartment_name,
                     hosp_compartment) {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity
    p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
    p[[gamma_name]] <- Inf
    p[[progression_name]] <- max(p[[progression_name]], 2)

    mod <- lancelot$new(p, 0, 1)

    info <- mod$info()
    state <- lancelot_initial(info, 1, p)

    # add individuals into the compartment
    if (hosp_compartment) {
      name_conf <- paste0(compartment_name, "_conf")
      name_unconf <- paste0(compartment_name, "_unconf")
      index_conf <- array(info$index[[name_conf]], info$dim[[name_conf]])
      index_unconf <- array(info$index[[name_unconf]], info$dim[[name_unconf]])
      state[index_conf] <- 50
      state[index_unconf] <- 50
    } else {
      index <- array(info$index[[compartment_name]],
                     info$dim[[compartment_name]])
      state[index] <- 50
    }

    mod$update_state(state = state)
    y <- mod$transform_variables(drop(mod$simulate(0:400)))

    if (hosp_compartment) {
      y[[compartment_name]] <- y[[name_conf]] + y[[name_unconf]]
    }

    z <- y[[compartment_name]]

    expect_true(any(z > 0))

    i <- seq_len(length(y$time) - 1L)
    if (length(dim(z)) == 5) {
      expect_equal(z[, , 2, , i + 1], z[, , 1, , i])
    } else {
      expect_equal(z[, , 2, i + 1], z[, , 1, i])
    }
  }

  helper("gamma_E_step", "k_E", "E", FALSE)
  helper("gamma_A_step", "k_A", "I_A", FALSE)
  helper("gamma_P_step", "k_P", "I_P", FALSE)
  helper("gamma_C_1_step", "k_C_1", "I_C_1", FALSE)
  helper("gamma_C_2_step", "k_C_2", "I_C_2", FALSE)
  helper("gamma_ICU_pre_step", "k_ICU_pre", "ICU_pre", TRUE)
  helper("gamma_H_R_step", "k_H_R", "H_R", TRUE)
  helper("gamma_H_D_step", "k_H_D", "H_D", TRUE)
  helper("gamma_ICU_W_R_step", "k_ICU_W_R", "ICU_W_R", TRUE)
  helper("gamma_ICU_W_D_step", "k_ICU_W_D", "ICU_W_D", TRUE)
  helper("gamma_ICU_D_step", "k_ICU_D", "ICU_D", TRUE)
  helper("gamma_G_D_step", "k_G_D", "G_D", FALSE)
  helper("gamma_W_R_step", "k_W_R", "W_R", TRUE)
  helper("gamma_W_D_step", "k_W_D", "W_D", TRUE)
  helper("gamma_sero_pre_1", "k_sero_pre_1", "T_sero_pre_1", FALSE)
  helper("gamma_sero_pos_1", "k_sero_pos_1", "T_sero_pos_1", FALSE)
  helper("gamma_sero_pre_2", "k_sero_pre_2", "T_sero_pre_2", FALSE)
  helper("gamma_sero_pos_2", "k_sero_pos_2", "T_sero_pos_2", FALSE)
  helper("gamma_PCR_pre_step", "k_PCR_pre", "T_PCR_pre", FALSE)
  helper("gamma_PCR_pos_step", "k_PCR_pos", "T_PCR_pos", FALSE)
})


test_that("setting a gamma to 0 results in no progression", {
  helper <- function(gamma_name, progression_name, compartment_name,
                     hosp_compartment) {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity
    p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
    p[[gamma_name]] <- 0
    p[[progression_name]] <- max(p[[progression_name]], 2)

    mod <- lancelot$new(p, 0, 1)

    info <- mod$info()
    state <- lancelot_initial(info, 1, p)

    # add individuals into the compartment (only the first progression stage)
    if (hosp_compartment) {
      name_conf <- paste0(compartment_name, "_conf")
      name_unconf <- paste0(compartment_name, "_unconf")
      index_conf <- array(info$index[[name_conf]], info$dim[[name_conf]])
      index_unconf <- array(info$index[[name_unconf]], info$dim[[name_unconf]])
      if (length(dim(index_conf)) == 4) {
        state[index_conf[, , 1, ]] <- 50
        state[index_unconf[, , 1, ]] <- 50
      } else {
        state[index_conf[, , 1]] <- 50
        state[index_unconf[, , 1]] <- 50
      }
    } else {
      index <- array(info$index[[compartment_name]],
                     info$dim[[compartment_name]])
      if (length(dim(index)) == 4) {
        state[index[, , 1, ]] <- 50
      } else {
        state[index[, , 1]] <- 50
      }
    }

    mod$update_state(state = state)
    y <- mod$transform_variables(drop(mod$simulate(0:400)))

    if (hosp_compartment) {
      y[[compartment_name]] <- y[[name_conf]] + y[[name_unconf]]
    }

    z <- y[[compartment_name]]

    expect_true(any(z > 0))

    if (length(dim(z)) == 5) {
      expect_true(all(z[, , 2, , ] == 0))
    } else {
      expect_true(all(z[, , 2, ] == 0))
    }
  }

  p <- lancelot_parameters(0, "england")
  helper("gamma_E_step", "k_E", "E", FALSE)
  helper("gamma_A_step", "k_A", "I_A", FALSE)
  helper("gamma_P_step", "k_P", "I_P", FALSE)
  helper("gamma_C_1_step", "k_C_1", "I_C_1", FALSE)
  helper("gamma_C_2_step", "k_C_2", "I_C_2", FALSE)
  helper("gamma_ICU_pre_step", "k_ICU_pre", "ICU_pre", TRUE)
  helper("gamma_H_R_step", "k_H_R", "H_R", TRUE)
  helper("gamma_H_D_step", "k_H_D", "H_D", TRUE)
  helper("gamma_ICU_W_R_step", "k_ICU_W_R", "ICU_W_R", TRUE)
  helper("gamma_ICU_W_D_step", "k_ICU_W_D", "ICU_W_D", TRUE)
  helper("gamma_ICU_D_step", "k_ICU_D", "ICU_D", TRUE)
  helper("gamma_G_D_step", "k_G_D", "G_D", FALSE)
  helper("gamma_W_R_step", "k_W_R", "W_R", TRUE)
  helper("gamma_W_D_step", "k_W_D", "W_D", TRUE)
  helper("gamma_sero_pre_1", "k_sero_pre_1", "T_sero_pre_1", FALSE)
  helper("gamma_sero_pos_1", "k_sero_pos_1", "T_sero_pos_1", FALSE)
  helper("gamma_sero_pre_2", "k_sero_pre_2", "T_sero_pre_2", FALSE)
  helper("gamma_sero_pos_2", "k_sero_pos_2", "T_sero_pos_2", FALSE)
  helper("gamma_PCR_pre_step", "k_PCR_pre", "T_PCR_pre", FALSE)
  helper("gamma_PCR_pos_step", "k_PCR_pos", "T_PCR_pos", FALSE)
})


test_that("No one is unconfirmed, if p_star = 1", {
  set.seed(1)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_star_step[, ] <- 1

  p$gamma_ICU_pre_step <- Inf
  p$gamma_H_R_step <- Inf
  p$gamma_H_D_step <- Inf
  p$gamma_ICU_W_R_step <- Inf
  p$gamma_ICU_W_D_step <- Inf
  p$gamma_ICU_D_step <- Inf
  p$gamma_W_R_step <- Inf
  p$gamma_W_D_step <- Inf

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p), time = 0)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(all(y$H_R_unconf == 0))
  expect_true(any(y$H_R_conf > 0))
  expect_true(all(y$H_D_unconf == 0))
  expect_true(any(y$H_D_conf > 0))
  expect_true(all(y$ICU_pre_unconf == 0))
  expect_true(any(y$ICU_pre_conf > 0))
  expect_true(all(y$ICU_W_R_unconf == 0))
  expect_true(any(y$ICU_W_R_conf > 0))
  expect_true(all(y$ICU_W_D_unconf == 0))
  expect_true(any(y$ICU_W_D_conf > 0))
  expect_true(all(y$ICU_D_unconf == 0))
  expect_true(any(y$ICU_D_conf > 0))
  expect_true(all(y$W_R_unconf == 0))
  expect_true(any(y$W_R_conf > 0))
  expect_true(all(y$W_D_unconf == 0))
  expect_true(any(y$W_D_conf > 0))

  admit_conf <- apply(y$H_R_conf[, 1, , , ] +
                        y$H_D_conf[, 1, , , ] +
                        y$ICU_pre_conf[, 1, , , ], 1, sum)

  expect_true(all(diff(y$cum_admit_conf) == admit_conf[-1]))
  expect_true(all(y$cum_new_conf == 0))
})


test_that("No one is confirmed, if p_star = 0 and gamma_U = 0", {
  ## RGF: I cannot replicate this failure on either of the Mac systems
  ## I have access to.
  skip_on_mac_gha()
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_star_step[, ] <- 0
  p$gamma_U_step <- 0

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  mod$update_state(state = lancelot_initial(info, 1, p), time = 0)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$H_R_unconf > 0))
  expect_true(all(y$H_R_conf == 0))
  expect_true(any(y$H_D_unconf > 0))
  expect_true(all(y$H_D_conf == 0))
  expect_true(any(y$ICU_pre_unconf > 0))
  expect_true(all(y$ICU_pre_conf == 0))
  expect_true(any(y$ICU_W_R_unconf > 0))
  expect_true(all(y$ICU_W_R_conf == 0))
  expect_true(any(y$ICU_W_D_unconf > 0))
  expect_true(all(y$ICU_W_D_conf == 0))
  expect_true(any(y$ICU_D_unconf > 0))
  expect_true(all(y$ICU_D_conf == 0))
  expect_true(any(y$W_R_unconf > 0))
  expect_true(all(y$W_R_conf == 0))
  expect_true(any(y$W_D_unconf > 0))
  expect_true(all(y$W_D_conf == 0))
  expect_true(all(y$admit_new_conf == 0))
  expect_true(all(y$cum_new_conf == 0))
})


test_that("Instant confirmation if p_star = 0 and gamma_U = Inf", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  p$p_star_step[, ] <- 0

  p$gamma_U_step <- Inf
  p$gamma_ICU_pre_step <- Inf
  p$gamma_H_R_step <- Inf
  p$gamma_H_D_step <- Inf
  p$gamma_ICU_W_R_step <- Inf
  p$gamma_ICU_W_D_step <- Inf
  p$gamma_ICU_D_step <- Inf
  p$gamma_W_R_step <- Inf
  p$gamma_W_D_step <- Inf

  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)

  ## We want to set ICU_W_R_unconf[, 1, ], ICU_W_D_unconf[, 1, ] and
  ## ICU_D_unconf[, 1, ] to 50
  y0[info$index$ICU_W_R_unconf[1:19]] <- 50
  y0[info$index$ICU_W_D_unconf[1:19]] <- 50
  y0[info$index$ICU_D_unconf[1:19]] <- 50

  mod$update_state(state = y0, time = 0)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))
  n <- length(y$time)

  ## Check hosp_R
  expect_true(all(y$H_R_conf[, , 1, , ] == 0))
  expect_equal(y$H_R_conf[, , 2, , -1], y$H_R_unconf[, , 1, , -n])
  expect_true(all(y$H_R_unconf[, , 2, , ] == 0))

  ## Check hosp_D
  expect_true(all(y$H_D_conf[, , 1, , ] == 0))
  expect_equal(y$H_D_conf[, , 2, , -1], y$H_D_unconf[, , 1, , -n])
  expect_true(all(y$H_D_unconf[, , 2, , ] == 0))

  ## Check ICU_pre
  expect_true(all(y$ICU_pre_conf[, , 1, , ] == 0))
  expect_equal(y$ICU_pre_conf[, , 2, , -1], y$ICU_pre_unconf[, , 1, , -n])
  expect_true(all(y$ICU_pre_unconf[, , 2, , ] == 0))

  ## Check ICU_D/ICU_S_R/ICU_S_D
  expect_equal(y$ICU_D_conf[, , 2, , 2], y$ICU_D_unconf[, , 1, , 1])
  expect_true(all(y$ICU_D_unconf[, , 2, , ] == 0))
  expect_equal(y$ICU_W_R_conf[, , 2, , 2], y$ICU_W_R_unconf[, , 1, , 1])
  expect_true(all(y$ICU_W_R_unconf[, , 2, , ] == 0))
  expect_true(all(y$ICU_W_D_unconf[, , 2, , ] == 0))
  expect_equal(y$ICU_W_D_conf[, , 2, , 2], y$ICU_W_D_unconf[, , 1, , 1])
  expect_equal(y$ICU_D_conf[, , 1, , -1] + y$ICU_W_R_conf[, , 1, , -1] +
                 y$ICU_W_D_conf[, , 1, , -1], y$ICU_pre_conf[, , 2, , -n])

  ## Check stepdown_R
  expect_equal(y$W_R_conf[, , 2, , 2],
               y$W_R_unconf[, , 1, , 1])
  expect_equal(y$W_R_conf[, , 1, , -1],
               y$ICU_W_R_conf[, , 2, , -n])
  expect_true(all(y$W_R_unconf[, , 2, , ] == 0))

  ## Check stepdown_D
  expect_equal(y$W_D_conf[, , 2, , 2],
               y$W_D_unconf[, , 1, , 1])
  expect_equal(y$W_D_conf[, , 1, , -1], y$ICU_W_D_conf[, , 2, , -n])
  expect_true(all(y$W_D_unconf[, , 2, , ] == 0))

  new_conf <- apply(y$H_R_conf[, , 2, , ] +
                      y$H_D_conf[, , 2, , ] +
                      y$ICU_pre_conf[, , 2, , ], 2, sum)
  new_conf[2] <- new_conf[2] +
    sum(y$ICU_W_R_conf[, , 2, , 2] +
          y$ICU_W_D_conf[, , 2, , 2] +
          y$ICU_D_conf[, , 2, , 2] +
          y$R_stepdown_conf[, , 2, , 2] +
          y$W_R_conf[, , 2, , 2])
  expect_true(all(diff(y$cum_new_conf) == new_conf[-1]))

  expect_true(all(y$cum_admit_conf == 0))
})


test_that("tots all summed correctly ", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))
  expect_true(all(y$general_tot == apply(y$ICU_pre_conf, 5, sum) +
                    apply(y$H_R_conf, 5, sum) +
                    apply(y$H_D_conf, 5, sum) +
                    apply(y$W_R_conf, 5, sum) +
                    apply(y$W_D_conf, 5, sum)))
  expect_true(all(y$ICU_tot == apply(y$ICU_W_R_conf, 5, sum) +
                    apply(y$ICU_W_D_conf, 5, sum) +
                    apply(y$ICU_D_conf, 5, sum)))
  expect_true(all(y$hosp_tot == y$ICU_tot + y$general_tot))
  expect_true(all(y$D_hosp_tot == apply(y$D_hosp, c(2, 3), sum)))
  expect_true(all(y$D_comm_tot == apply(y$D_non_hosp[1:18, ], 2, sum)))
  expect_true(all(y$D_lancelot_tot == y$D_non_hosp[19, ]))
  expect_true(all(y$D_tot == y$D_hosp_tot + y$D_lancelot_tot + y$D_comm_tot))

  # check the positivity sums
  expect_true(all(y$sero_pos_1 == apply(y$T_sero_pos_1[4:13, , , 1, ], 3, sum)))
  expect_true(all(y$sero_pos_2 == apply(y$T_sero_pos_2[4:13, , , 1, ], 3, sum)))
  expect_true(all(y$react_pos == apply(y$T_PCR_pos[2:18, , , 1, ], 3, sum)))
  expect_true(all(y$ons_pos == apply(y$T_PCR_pos[1, , , 1, ], 2, sum) * 3 / 5 +
                    apply(y$T_PCR_pos[2:18, , , 1, ], 3, sum)))
})


test_that("Symptomatic cases by age add up correctly", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$sympt_cases_inc ==
                    y$sympt_cases_under15_inc + y$sympt_cases_15_24_inc +
                    y$sympt_cases_25_49_inc + y$sympt_cases_50_64_inc +
                    y$sympt_cases_65_79_inc + y$sympt_cases_80_plus_inc))

  expect_true(all(y$sympt_cases_over25_inc ==
                    y$sympt_cases_25_49_inc + y$sympt_cases_50_64_inc +
                    y$sympt_cases_65_79_inc + y$sympt_cases_80_plus_inc))
})


test_that("Infections and hospitlisations incidence by age add up correctly", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$infections_inc == colSums(y$infections_inc_age)))

  expect_true(all(y$hospitalisations_inc ==
                    colSums(y$hospitalisations_inc_age)))
})


test_that("Disaggregated and aggregated data streams add up correctly", {
  p <- lancelot_parameters(0, "england")
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(round(y$D_hosp_inc) ==
                    round(y$D_hosp_0_49_inc +
                            y$D_hosp_50_54_inc + y$D_hosp_55_59_inc +
                            y$D_hosp_60_64_inc + y$D_hosp_65_69_inc +
                            y$D_hosp_70_74_inc + y$D_hosp_75_79_inc +
                            y$D_hosp_80_plus_inc)))
  expect_true(all(round(y$D_comm_inc) ==
                    round(y$D_comm_0_49_inc +
                            y$D_comm_50_54_inc + y$D_comm_55_59_inc +
                            y$D_comm_60_64_inc + y$D_comm_65_69_inc +
                            y$D_comm_70_74_inc + y$D_comm_75_79_inc +
                            y$D_comm_80_plus_inc)))
  expect_true(all(round(y$admit_conf_inc + y$new_conf_inc) ==
                    round(y$all_admission_0_9_conf_inc +
                            y$all_admission_10_19_conf_inc +
                            y$all_admission_20_29_conf_inc +
                            y$all_admission_30_39_conf_inc +
                            y$all_admission_40_49_conf_inc +
                            y$all_admission_50_59_conf_inc +
                            y$all_admission_60_69_conf_inc +
                            y$all_admission_70_79_conf_inc +
                            y$all_admission_80_plus_conf_inc)))
  expect_true(all(round(y$react_pos) ==
                    round(y$react_5_24_pos +
                            y$react_25_34_pos +
                            y$react_35_44_pos +
                            y$react_45_54_pos +
                            y$react_55_64_pos +
                            y$react_65_plus_pos)))
})


test_that("Individuals cannot infect in compartment with zero transmission", {
  helper <- function(transmission_name, compartment_name, gamma_name) {
    ## Use a large beta so that infections would be immediate
    ## No seeding
    p <- lancelot_parameters(0, "england", beta_value = 1e9,
                             initial_seed_size = 0)

    ## set all transmission parameters to 1
    p$I_A_transmission <- 1
    p$I_P_transmission <- 1
    p$I_C_1_transmission <- 1
    p$I_C_2_transmission <- 1
    p$hosp_transmission <- 1
    p$ICU_transmission <- 1
    p$G_D_transmission <- 1

    ## set transmission parameters being tested to 0
    p[[transmission_name]] <- 0

    ## set relevant gamma to 0 so no progression
    p[[gamma_name]] <- 0

    mod <- lancelot$new(p, 0, 1)

    info <- mod$info()
    y0 <- lancelot_initial(info, 1, p)

    ## put individuals in compartment being tested
    index <- info$index[[compartment_name]]
    y0[index] <- 50

    mod$update_state(state = y0)
    y <- mod$transform_variables(
      drop(mod$simulate(seq(0, 400, by = 4))))

    ## Susceptible population is never drawn down:
    expect_equal(y$S, array(y$S[, , 1], c(nrow(y$S), 1, 101)))
    ## Ignore group 4 where we weight 1 if all are 0
    expect_true(all(y$I_weighted[-4, , , ] == 0))
  }

  helper("I_A_transmission", "I_A", "gamma_A_step")
  helper("I_P_transmission", "I_P", "gamma_P_step")
  helper("I_C_1_transmission", "I_C_1", "gamma_C_1_step")
  helper("I_C_2_transmission", "I_C_2", "gamma_C_2_step")
  helper("G_D_transmission", "G_D", "gamma_G_D_step")
  helper("hosp_transmission", "H_D_unconf", "gamma_H_D_step")
  helper("hosp_transmission", "H_D_conf", "gamma_H_D_step")
  helper("hosp_transmission", "H_R_unconf", "gamma_H_R_step")
  helper("hosp_transmission", "H_R_conf", "gamma_H_R_step")
  helper("hosp_transmission", "ICU_pre_unconf", "gamma_ICU_pre_step")
  helper("hosp_transmission", "ICU_pre_conf", "gamma_ICU_pre_step")
  helper("ICU_transmission", "ICU_D_unconf", "gamma_ICU_D_step")
  helper("ICU_transmission", "ICU_D_conf", "gamma_ICU_D_step")
  helper("ICU_transmission", "ICU_W_D_unconf", "gamma_ICU_W_D_step")
  helper("ICU_transmission", "ICU_W_D_conf", "gamma_ICU_W_D_step")
  helper("ICU_transmission", "ICU_W_R_unconf", "gamma_ICU_W_R_step")
  helper("ICU_transmission", "ICU_W_R_conf", "gamma_ICU_W_R_step")
})


test_that("Individuals can infect in compartment with non-zero transmission", {
  helper <- function(transmission_name, compartment_name, gamma_name) {
    ## Use a large beta so that infections would be immediate
    ## No seeding
    p <- lancelot_parameters(0, "england", beta_value = 1e9,
                             initial_seed_size = 0)

    ## set all transmission parameters to 0
    p$I_A_transmission <- 0
    p$I_P_transmission <- 0
    p$I_C_1_transmission <- 0
    p$I_C_2_transmission <- 0
    p$hosp_transmission <- 0
    p$ICU_transmission <- 0
    p$G_D_transmission <- 0

    ## set transmission parameter being tested to a non-zero value
    p[[transmission_name]] <- 0.9

    ## set relevant gamma to 0 so no progression
    p[[gamma_name]] <- 0

    mod <- lancelot$new(p, 0, 1)

    info <- mod$info()
    y0 <- lancelot_initial(info, 1, p)

    ## put individuals in compartment being tested
    index <- info$index[[compartment_name]]
    y0[index] <- 50
    index_I_weighted <- info$index$I_weighted
    y0[index_I_weighted] <- p[[transmission_name]] * 50 *
      info$dim[[compartment_name]][[3]]

    mod$update_state(state = y0)
    y <- mod$transform_variables(
      drop(mod$simulate(c(0, 1))))

    ## Susceptible population is immediately infected:
    expect_true(all(y$S[, , 2] == 0))
    ## I_weighted calculated as expected
    expect_true(all(y$I_weighted[, , , 2] == y$I_weighted[, , , 1]))
  }

  helper("I_A_transmission", "I_A", "gamma_A_step")
  helper("I_P_transmission", "I_P", "gamma_P_step")
  helper("I_C_1_transmission", "I_C_1", "gamma_C_1_step")
  helper("I_C_2_transmission", "I_C_2", "gamma_C_2_step")
  helper("G_D_transmission", "G_D", "gamma_G_D_step")
  helper("hosp_transmission", "H_D_unconf", "gamma_H_D_step")
  helper("hosp_transmission", "H_D_conf", "gamma_H_D_step")
  helper("hosp_transmission", "H_R_unconf", "gamma_H_R_step")
  helper("hosp_transmission", "H_R_conf", "gamma_H_R_step")
  helper("hosp_transmission", "ICU_pre_unconf", "gamma_ICU_pre_step")
  helper("hosp_transmission", "ICU_pre_conf", "gamma_ICU_pre_step")
  helper("ICU_transmission", "ICU_D_unconf", "gamma_ICU_D_step")
  helper("ICU_transmission", "ICU_D_conf", "gamma_ICU_D_step")
  helper("ICU_transmission", "ICU_W_D_unconf", "gamma_ICU_W_D_step")
  helper("ICU_transmission", "ICU_W_D_conf", "gamma_ICU_W_D_step")
  helper("ICU_transmission", "ICU_W_R_unconf", "gamma_ICU_W_R_step")
  helper("ICU_transmission", "ICU_W_R_conf", "gamma_ICU_W_R_step")
})


test_that("Severity probabilities correctly capped at 1", {
  ## We're going to try a number of very similar tests here, so a
  ## helper function will help run a model with a severity probability set to
  ## 1 and then also set to 1.5. Due to capping at 1 in the model, we
  ## expect these to produce the same results with the same seed
  helper <- function(prob_name) {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity

    p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
    ## First set value to 1
    p[[paste0(prob_name, "_step")]][, ] <- 1

    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    mod$update_state(state = lancelot_initial(info, 1, p))
    y <- mod$simulate(seq(0, 400, by = 4))

    ## Now set value to 1.5
    p[[paste0(prob_name, "_step")]][, ] <- 1.5

    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    mod$update_state(state = lancelot_initial(info, 1, p))
    y2 <- mod$simulate(seq(0, 400, by = 4))

    expect_true(all(y == y2, na.rm = TRUE))
  }

  helper("p_C")
  helper("p_H")
  helper("p_ICU")
  helper("p_ICU_D")
  helper("p_H_D")
  helper("p_W_D")
  helper("p_G_D")
  helper("p_R")

})


test_that("Severity outputs are correctly calculated", {

  helper <- function(prob_C, prob_H, prob_G_D,
                     prob_ICU, prob_ICU_D, prob_H_D, prob_W_D,
                     ihr, ifr, hfr,
                     ihr_strain, ifr_strain, hfr_strain) {

    p <- lancelot_parameters(0, "england", waning_rate = 1 / 20)
    p$p_C_step[, ] <- prob_C %||% p$p_C_step
    p$p_H_step[, ] <- prob_H %||% p$p_H_step
    p$p_G_D_step[, ] <- prob_G_D %||% p$p_G_D_step
    p$p_ICU_step[, ] <- prob_ICU %||% p$p_ICU_step
    p$p_ICU_D_step[, ] <- prob_ICU_D %||% p$p_ICU_D_step
    p$p_H_D_step[, ] <- prob_H_D %||% p$p_H_D_step
    p$p_W_D_step[, ] <- prob_W_D %||% p$p_W_D_step

    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    mod$update_state(state = lancelot_initial(info, 1, p))
    y <- mod$transform_variables(
      drop(mod$simulate(seq(0, 400, by = 4))))

    mod_ifr <- y$ifr[!is.na(y$ifr)][-1]
    mod_ihr <- y$ihr[!is.na(y$ihr)][-1]
    mod_hfr <- y$hfr[!is.na(y$hfr)][-1]

    mod_ihr_strain <- y$ihr_strain[!is.na(y$ihr_strain)][-1]
    mod_ifr_strain <- y$ifr_strain[!is.na(y$ifr_strain)][-1]
    mod_hfr_strain <- y$hfr_strain[!is.na(y$hfr_strain)][-1]

    if (length(mod_ifr) > 0) {
      expect_vector_equal(mod_ifr, ifr)
    }
    if (length(mod_ihr) > 0) {
      expect_vector_equal(mod_ihr, ihr)
    }
    if (length(mod_hfr) > 0) {
      expect_vector_equal(mod_hfr, hfr)
    }
    if (length(mod_ifr_strain) > 0) {
      expect_vector_equal(mod_ifr_strain, ifr_strain)
    }
    if (length(mod_ihr_strain) > 0) {
      expect_vector_equal(mod_ihr_strain, ihr_strain)
    }
    if (length(mod_hfr_strain) > 0) {
      expect_vector_equal(mod_hfr_strain, hfr_strain)
    }
  }

  # Test p_C = 0
  helper(0, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0)

  # Test p_C = 1 & p_H = 0
  helper(1, 0, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0)

  # Test p_C = 1 & p_H = 1 & p_G_D = 1
  helper(1, 1, 1, NULL, NULL, NULL, NULL, 0, 1, 0, 0, 1, 0)

  # Test p_C = 1 & p_H = 1 & p_G_D = 0 & p_ICU = 0 & p_H_D = 0
  helper(1, 1, 0, 0, NULL, 0, NULL, 1, 0, 0, 1, 0, 0)

  # Test p_C = 1 & p_H = 1 & p_G_D = 0 & p_ICU = 0 & p_H_D = 1
  helper(1, 1, 0, 0, NULL, 1, NULL, 1, 1, 1, 1, 1, 1)

  # Test p_C = 1 & p_H = 1 & p_G_D = 0 & p_ICU = 1 & p_ICU_D = 0 & p_W_D = 1
  helper(1, 1, 0, 1, 0, NULL, 1, 1, 1, 1, 1, 1, 1)

  # Test p_C = 1 & p_H = 1 & p_G_D = 0 & p_ICU = 1 & p_ICU_D = 0 & p_W_D = 0
  helper(1, 1, 0, 1, 0, NULL, 0, 1, 0, 0, 1, 0, 0)

  # Test p_C = 1 & p_H = 1 & p_G_D = 0 & p_ICU = 1 & p_ICU_D = 1
  helper(1, 1, 0, 1, 1, NULL, NULL, 1, 1, 1, 1, 1, 1)
})


test_that("Severity by age is calculated parametrically", {

  helper <- function(p, i) {
    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    state <- lancelot_initial(info, 1, p)

    index_i <- array(info$index[[i]], info$dim[[i]])

    mod$update_state(state = state)
    mod$set_index(info$index[[i]])

    y <- mod$simulate(seq(0, 800, by = 4))

    expect_equal(length(y), prod(info$dim[[i]]) * 201)

    y <- array(y, c(info$dim[[i]], 201))
  }

  p <- lancelot_parameters(1, "uk", carehome_beds = 0)

  ihr <- p$p_C_step * p$p_H_step * (1 - p$p_G_D_step)
  hfr <- (1 - p$p_ICU_step) * p$p_H_D_step +
    p$p_ICU_step * (p$p_ICU_D_step + (1 - p$p_ICU_D_step) * p$p_W_D_step)
  ifr <- ihr * hfr + p$p_C_step * p$p_H_step * p$p_G_D_step

  # There are extremely low discrepancies due to integer precision, in random
  # number generation. We'll round to 10 digits
  y <- helper(p, "ifr_age")[c(1:17), 201]
  x <- ifr[which(!is.na(y))]
  expect_vector_equal(x, y, 10)

  y <- helper(p, "ihr_age")[c(1:17), 101]
  x <- ihr[which(!is.na(y))]
  expect_vector_equal(x, y, 10)

  # We have some NAs in young age bands in y - we'll ignore
  y <- helper(p, "hfr_age")[c(1:17), 101]
  x <- hfr[which(!is.na(y))]
  y <- y[which(!is.na(y))]
  expect_vector_equal(x, y, 10)
})


test_that("Effective susceptible and protected calculation work as expected", {
  p <- lancelot_parameters(0, "uk")
  mod <- lancelot$new(p, 0, 1)
  info <- mod$info()
  y0 <- lancelot_initial(info, 1, p)
  mod$update_state(state = lancelot_initial(info, 1, p))
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## With single strain and no vaccination, effective_susceptible should just
  ## equal the total susceptible
  expect_true(all(apply(y$S, 3, sum) == y$effective_susceptible))
  ## All individuals in R should be fully protected
  expect_true(all(apply(y$R, 4, sum) == y$protected_R_unvaccinated))
  ## No-one should be vaccinated
  expect_true(all(y$protected_S_vaccinated == 0))
  expect_true(all(y$protected_R_vaccinated == 0))
})

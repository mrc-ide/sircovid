context("carehomes (check)")

test_that("N_tot, N_tot2 and N_tot3 stay constant", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "uk", waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  ## This is not quite correct, and I don't really know why. I think
  ## that this is the rounding that we're doing to shuffle the
  ## carehome residents around. However, it does mean that this is
  ## different between runs.  Do a careful check with the sircovid
  ## model. We are off by one! individual.
  expect_true(all(y$N_tot3 - mod$transform_variables(y0)$N_tot3 == 0))
  expect_true(all(y$N_tot2 - mod$transform_variables(y0)$N_tot2 == 0))
  expect_true(all(y$N_tot - mod$transform_variables(y0)$N_tot == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot2 == 0))
  expect_true(all(colSums(y$N_tot) - y$N_tot3 == 0))
})


test_that("there are no infections when beta is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  waning_rate <- rep(1 / 20, 19)
  waning_rate[4] <- 0 # no waning in group with seeded infections
  # otherwise S can go up as these infected individuals loose immunity

  p <- carehomes_parameters(0, "england", beta_value = 0,
                            waning_rate = waning_rate)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(nrow(s), 1, 101)))
})


test_that("everyone is infected when beta is large", {
  p <- carehomes_parameters(0, "england", beta_value = 1e9)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(
    mod$simulate(seq(0, 400, by = 4))))
  expect_true(all(y$S[, 1, -1] == 0))
})


test_that("noone stays in R, T_sero_neg or T_PCR_neg if waning rate is very
          large", {
  # with a large waning rate and beta = 0,
  # people can move from R to S but not outside of S
  # therefore R should quickly get empty (and T_sero_neg and T_PCR_neg as well)
  p <- carehomes_parameters(0, "england",
                            beta_value = 0, # to forbid movement out of S
                            waning_rate = Inf) # to force movement out of R
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  state <- carehomes_initial(info, 1, p)$state

  # move everyone to R, T_sero_neg and T_PCR_neg initially
  index_S <- array(info$index$S, info$dim$S)
  index_R <- array(info$index$R, info$dim$R)
  index_T_sero_neg <- array(info$index$T_sero_neg, info$dim$T_sero_neg)
  index_T_PCR_neg <- array(info$index$T_PCR_neg, info$dim$T_PCR_neg)
  state[index_R] <- rowSums(array(state[index_S], info$dim$S))
  state[index_T_sero_neg] <- rowSums(array(state[index_S], info$dim$S))
  state[index_T_PCR_neg] <- rowSums(array(state[index_S], info$dim$S))
  state[index_S] <- 0

  mod$set_state(state)
  y <- mod$transform_variables(drop(
    mod$simulate(seq(0, 400, by = 4))))

  # other than in the 4th age group (where infections are seeded)
  # after the first day (4 times steps), R is empty
  expect_true(all(y$R[-4, , -1, ] == 0))
  # so is T_sero_neg
  expect_true(all(y$T_sero_neg[-4, , -1, ] == 0))
  # so is T_PCR_neg
  expect_true(all(y$T_PCR_neg[-4, , -1, ] == 0))

})

test_that("R, T_sero_neg and T_PCR_neg are all non-decreasing and S is
          non-increasing if waning rate is 0", {
  p <- carehomes_parameters(0, "england",
                            waning_rate = 0)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(
    mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(diff(t(drop(y$R))) >= 0))
  expect_true(all(diff(t(drop(y$T_sero_neg))) >= 0))
  expect_true(all(diff(t(drop(y$T_PCR_neg))) >= 0))
  expect_true(all(diff(t(y$S[, 1, ])) <= 0))

})


test_that("No one is infected if I and E are 0 at t = 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y <- carehomes_initial(info, 1, p)$state
  y[info$index$I_A] <- 0
  y[info$index$T_sero_pre] <- 0
  y[info$index$T_PCR_pos] <- 0

  mod$set_state(y)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(nrow(s), 1, 101)))
})


test_that("No one is hospitalised, no-one dies if p_C is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$p_C[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
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


test_that("No one is hospitalised, no-one dies if psi_H is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$psi_H[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
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
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$p_H_step <- 1
  p$psi_H[] <- 1
  p$p_G_D_step <- 1
  p$psi_G_D[] <- 1
  p$p_C[] <- 1

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


test_that("No one is hospitalised, no-one recovers in edge case 2", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$p_H_step <- 1
  p$psi_H[] <- 1
  p$p_G_D_step <- 1
  p$psi_G_D[] <- 1
  p$p_C[] <- 1

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


test_that("No one dies in the community if psi_G_D is 0", {
  set.seed(1)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$psi_G_D[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_P > 0))
  expect_true(any(y$I_C_1 > 0))
  expect_true(any(y$I_C_2 > 0))
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
    p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
    p$p_ICU_step <- ifelse(is.null(prob_ICU),
                                    p$p_ICU_step, 1)
    p$psi_ICU[] <- prob_ICU %||% p$psi_ICU[]
    p$p_H_D_step <- ifelse(is.null(prob_H_D),
                                 p$p_H_D_step, 1)
    p$psi_H_D[] <- prob_H_D %||% p$psi_H_D[]
    p$p_ICU_D_step <- ifelse(is.null(prob_ICU_D),
                                 p$p_ICU_D_step, 1)
    p$psi_ICU_D[] <- prob_ICU_D %||% p$psi_ICU_D[]
    p$p_W_D_step <- ifelse(is.null(prob_W_D),
                                 p$p_W_D_step, 1)
    p$psi_W_D[] <- prob_W_D %||% p$psi_W_D[]

    mod <- carehomes$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    mod$set_state(carehomes_initial(info, 1, p)$state)
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
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$p_sero_pos[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$T_sero_neg > 0))
  expect_true(all(y$T_sero_pos == 0))
})


test_that("No one does not seroconvert and no one seroreverts
          if p_sero_pos is 1 and gamma_sero_pos is 0", {
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)

  p$p_sero_pos[] <- 1
  ## set gamma_sero_pos to 0 so no-one seroreverts
  p$gamma_sero_pos <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(all(y$T_sero_neg == 0))
  expect_true(any(y$T_sero_pos > 0))
})


test_that("T_sero_pre parameters work as expected", {
  helper <- function(p_sero_pre_1, gamma_sero_pre_1, gamma_sero_pre_2) {
    p <- carehomes_parameters(0, "uk")
    p$p_sero_pre_1 <- p_sero_pre_1
    p$gamma_sero_pre_1 <- gamma_sero_pre_1
    p$gamma_sero_pre_2 <- gamma_sero_pre_2

    mod <- carehomes$new(p, 0, 1)
    info <- mod$info()

    y0 <- carehomes_initial(info, 1, p)$state
    y0[info$index$T_sero_pre] <- 0

    mod$set_state(y0)
    mod$transform_variables(drop(mod$simulate(0:400)))
  }

  ## p_sero_pre_1 = 1, expect no cases in T_sero_pre_2 stream
  y <- helper(1, 1, 0.5)
  expect_true(all(y$T_sero_pre[, , 2, , ] == 0))

  ## p_sero_pre_1 = 0, expect no cases in T_sero_pre_1 stream
  y <- helper(0, 1, 0.5)
  expect_true(all(y$T_sero_pre[, , 1, , ] == 0))

  ## gamma_sero_pre_1 = gamma_sero_pre_2 = 0, expect no cases in T_sero_pos
  y <- helper(0.5, 0, 0)
  expect_true(all(y$T_sero_pos == 0))
  expect_true(all(y$T_sero_neg == 0))

  ## gamma_sero_pre_1 = Inf, gamma_sero_pre_2 = 0, expect progression in one
  ## time-step to T_sero_neg/T_sero_pos just from T_sero_pre_1
  y <- helper(0.5, Inf, 0)
  n <- length(y$time)
  expect_equal(diff(t(apply(y$T_sero_pos, c(1, 5), sum) + drop(y$T_sero_neg))),
               t(y$T_sero_pre[, , 1, 1, -n]))

  ## gamma_sero_pre_1 = 0, gamma_sero_pre_2 = Inf, expect progression in one
  ## time-step to T_sero_neg/T_sero_pos just from T_sero_pre_2
  y <- helper(0.5, 0, Inf)
  n <- length(y$time)
  expect_equal(diff(t(apply(y$T_sero_pos, c(1, 5), sum) + drop(y$T_sero_neg))),
               t(y$T_sero_pre[, , 2, 1, -n]))

  ## gamma_sero_pre_1 = Inf, gamma_sero_pre_2 = Inf, expect progression in one
  ## time-step to T_sero_neg/T_sero_pos from both T_sero_pre_1 and T_sero_pre_2
  y <- helper(0.5, Inf, Inf)
  n <- length(y$time)
  expect_equal(diff(t(apply(y$T_sero_pos, c(1, 5), sum) + drop(y$T_sero_neg))),
               t(apply(y$T_sero_pre[, , , 1, -n], c(1, 3), sum)))
})


test_that("setting a gamma to Inf results immediate progression", {
  helper <- function(gamma_name, progression_name, compartment_name,
                     hosp_compartment) {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity
    p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
    p[[gamma_name]] <- Inf
    p[[progression_name]] <- max(p[[progression_name]], 2)

    mod <- carehomes$new(p, 0, 1)

    info <- mod$info()
    state <- carehomes_initial(info, 1, p)$state

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

    mod$set_state(state)
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

  helper("gamma_E", "k_E", "E", FALSE)
  helper("gamma_A", "k_A", "I_A", FALSE)
  helper("gamma_P", "k_P", "I_P", FALSE)
  helper("gamma_C_1", "k_C_1", "I_C_1", FALSE)
  helper("gamma_C_2", "k_C_2", "I_C_2", FALSE)
  helper("gamma_ICU_pre", "k_ICU_pre", "ICU_pre", TRUE)
  helper("gamma_H_R", "k_H_R", "H_R", TRUE)
  helper("gamma_H_D", "k_H_D", "H_D", TRUE)
  helper("gamma_ICU_W_R", "k_ICU_W_R", "ICU_W_R", TRUE)
  helper("gamma_ICU_W_D", "k_ICU_W_D", "ICU_W_D", TRUE)
  helper("gamma_ICU_D", "k_ICU_D", "ICU_D", TRUE)
  helper("gamma_G_D", "k_G_D", "G_D", FALSE)
  helper("gamma_W_R", "k_W_R", "W_R", TRUE)
  helper("gamma_W_D", "k_W_D", "W_D", TRUE)
  helper("gamma_sero_pos", "k_sero_pos", "T_sero_pos", FALSE)
  helper("gamma_PCR_pre", "k_PCR_pre", "T_PCR_pre", FALSE)
  helper("gamma_PCR_pos", "k_PCR_pos", "T_PCR_pos", FALSE)
})


test_that("setting a gamma to 0 results in no progression", {
  helper <- function(gamma_name, progression_name, compartment_name,
                     hosp_compartment) {
    ## waning_rate default is 0, setting to a non-zero value so that this test
    ## passes with waning immunity
    p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
    p[[gamma_name]] <- 0
    p[[progression_name]] <- max(p[[progression_name]], 2)

    mod <- carehomes$new(p, 0, 1)

    info <- mod$info()
    state <- carehomes_initial(info, 1, p)$state

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

    mod$set_state(state)
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

  p <- carehomes_parameters(0, "england")
  helper("gamma_E", "k_E", "E", FALSE)
  helper("gamma_A", "k_A", "I_A", FALSE)
  helper("gamma_P", "k_P", "I_P", FALSE)
  helper("gamma_C_1", "k_C_1", "I_C_1", FALSE)
  helper("gamma_C_2", "k_C_2", "I_C_2", FALSE)
  helper("gamma_ICU_pre", "k_ICU_pre", "ICU_pre", TRUE)
  helper("gamma_H_R", "k_H_R", "H_R", TRUE)
  helper("gamma_H_D", "k_H_D", "H_D", TRUE)
  helper("gamma_ICU_W_R", "k_ICU_W_R", "ICU_W_R", TRUE)
  helper("gamma_ICU_W_D", "k_ICU_W_D", "ICU_W_D", TRUE)
  helper("gamma_ICU_D", "k_ICU_D", "ICU_D", TRUE)
  helper("gamma_G_D", "k_G_D", "G_D", FALSE)
  helper("gamma_W_R", "k_W_R", "W_R", TRUE)
  helper("gamma_W_D", "k_W_D", "W_D", TRUE)
  helper("gamma_sero_pos", "k_sero_pos", "T_sero_pos", FALSE)
  helper("gamma_PCR_pre", "k_PCR_pre", "T_PCR_pre", FALSE)
  helper("gamma_PCR_pos", "k_PCR_pos", "T_PCR_pos", FALSE)
})


test_that("No one is unconfirmed, if p_star = 1", {
  set.seed(1)
  ## waning_rate default is 0, setting to a non-zero value so that this test
  ## passes with waning immunity
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$p_star_step <- 1
  p$gamma_H_R_step <- 1
  p$psi_star[] <- 1
  p$gamma_ICU_pre <- Inf
  p$gamma_H_R <- Inf
  p$gamma_H_D <- Inf
  p$gamma_ICU_W_R <- Inf
  p$gamma_ICU_W_D <- Inf
  p$gamma_ICU_D <- Inf
  p$gamma_W_R <- Inf
  p$gamma_W_D <- Inf

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state, 0)
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
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$psi_star[] <- 0
  p$gamma_U <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state, 0)
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
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  p$psi_star[] <- 0
  p$gamma_H_R_step <- 1
  p$gamma_U <- Inf
  p$gamma_ICU_pre <- Inf
  p$gamma_H_R <- Inf
  p$gamma_H_D <- Inf
  p$gamma_ICU_W_R <- Inf
  p$gamma_ICU_W_D <- Inf
  p$gamma_ICU_D <- Inf
  p$gamma_W_R <- Inf
  p$gamma_W_D <- Inf

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state

  ## We want to set ICU_W_R_unconf[, 1, ], ICU_W_D_unconf[, 1, ] and
  ## ICU_D_unconf[, 1, ] to 50
  y0[info$index$ICU_W_R_unconf[1:19]] <- 50
  y0[info$index$ICU_W_D_unconf[1:19]] <- 50
  y0[info$index$ICU_D_unconf[1:19]] <- 50

  mod$set_state(y0, 0)
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
  p <- carehomes_parameters(0, "england", waning_rate = 1 / 20)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
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
  expect_true(all(y$D_hosp_tot == apply(y$D_hosp, 2, sum)))
  expect_true(all(y$D_comm_tot == apply(y$D_non_hosp[1:18, ], 2, sum)))
  expect_true(all(y$D_carehomes_tot == y$D_non_hosp[19, ]))
  expect_true(all(y$D_tot == y$D_hosp_tot + y$D_carehomes_tot + y$D_comm_tot))

  # check the positivity sums
  expect_true(all(y$sero_pos == apply(y$T_sero_pos[4:13, , , 1, ], 3, sum)))
  expect_true(all(y$react_pos == apply(y$T_PCR_pos[2:18, , , 1, ], 3, sum)))
})


test_that("Individuals cannot infect in compartment with zero transmission", {
  helper <- function(transmission_name, compartment_name, gamma_name) {
    ## Use a large beta so that infections would be immediate
    p <- carehomes_parameters(0, "england", beta_value = 1e9)

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

    mod <- carehomes$new(p, 0, 1)

    info <- mod$info()
    y0 <- carehomes_initial(info, 1, p)$state

    ## remove initial asymptomatic individuals
    index_I_A <- info$index$I_A
    y0[index_I_A] <- 0
    index_I_weighted <- info$index$I_weighted
    y0[index_I_weighted] <- 0

    ## put individuals in compartment being tested
    index <- info$index[[compartment_name]]
    y0[index] <- 50

    mod$set_state(y0)
    y <- mod$transform_variables(
      drop(mod$simulate(seq(0, 400, by = 4))))

    ## Susceptible population is never drawn down:
    expect_equal(y$S, array(y$S[, , 1], c(nrow(y$S), 1, 101)))
    expect_true(all(y$I_weighted == 0))
  }

  helper("I_A_transmission", "I_A", "gamma_A")
  helper("I_P_transmission", "I_P", "gamma_P")
  helper("I_C_1_transmission", "I_C_1", "gamma_C_1")
  helper("I_C_2_transmission", "I_C_2", "gamma_C_2")
  helper("G_D_transmission", "G_D", "gamma_G_D")
  helper("hosp_transmission", "H_D_unconf", "gamma_H_D")
  helper("hosp_transmission", "H_D_conf", "gamma_H_D")
  helper("hosp_transmission", "H_R_unconf", "gamma_H_R")
  helper("hosp_transmission", "H_R_conf", "gamma_H_R")
  helper("hosp_transmission", "ICU_pre_unconf", "gamma_ICU_pre")
  helper("hosp_transmission", "ICU_pre_conf", "gamma_ICU_pre")
  helper("ICU_transmission", "ICU_D_unconf", "gamma_ICU_D")
  helper("ICU_transmission", "ICU_D_conf", "gamma_ICU_D")
  helper("ICU_transmission", "ICU_W_D_unconf", "gamma_ICU_W_D")
  helper("ICU_transmission", "ICU_W_D_conf", "gamma_ICU_W_D")
  helper("ICU_transmission", "ICU_W_R_unconf", "gamma_ICU_W_R")
  helper("ICU_transmission", "ICU_W_R_conf", "gamma_ICU_W_R")
})


test_that("Individuals can infect in compartment with non-zero transmission", {
  helper <- function(transmission_name, compartment_name, gamma_name) {
    ## Use a large beta so that infections would be immediate
    p <- carehomes_parameters(0, "england", beta_value = 1e9)

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

    mod <- carehomes$new(p, 0, 1)

    info <- mod$info()
    y0 <- carehomes_initial(info, 1, p)$state

    ## remove initial asymptomatic individuals
    index_I_A <- info$index$I_A
    y0[index_I_A] <- 0

    ## put individuals in compartment being tested
    index <- info$index[[compartment_name]]
    y0[index] <- 50
    index_I_weighted <- info$index$I_weighted
    y0[index_I_weighted] <- p[[transmission_name]] * 50 *
      info$dim[[compartment_name]][[3]]

    mod$set_state(y0)
    y <- mod$transform_variables(
      drop(mod$simulate(c(0, 1))))

    ## Susceptible population is immediately infected:
    expect_true(all(y$S[, , 2] == 0))
    ## I_weighted calculated as expected
    expect_true(all(y$I_weighted[, , 2] == y$I_weighted[, , 1]))
  }

  helper("I_A_transmission", "I_A", "gamma_A")
  helper("I_P_transmission", "I_P", "gamma_P")
  helper("I_C_1_transmission", "I_C_1", "gamma_C_1")
  helper("I_C_2_transmission", "I_C_2", "gamma_C_2")
  helper("G_D_transmission", "G_D", "gamma_G_D")
  helper("hosp_transmission", "H_D_unconf", "gamma_H_D")
  helper("hosp_transmission", "H_D_conf", "gamma_H_D")
  helper("hosp_transmission", "H_R_unconf", "gamma_H_R")
  helper("hosp_transmission", "H_R_conf", "gamma_H_R")
  helper("hosp_transmission", "ICU_pre_unconf", "gamma_ICU_pre")
  helper("hosp_transmission", "ICU_pre_conf", "gamma_ICU_pre")
  helper("ICU_transmission", "ICU_D_unconf", "gamma_ICU_D")
  helper("ICU_transmission", "ICU_D_conf", "gamma_ICU_D")
  helper("ICU_transmission", "ICU_W_D_unconf", "gamma_ICU_W_D")
  helper("ICU_transmission", "ICU_W_D_conf", "gamma_ICU_W_D")
  helper("ICU_transmission", "ICU_W_R_unconf", "gamma_ICU_W_R")
  helper("ICU_transmission", "ICU_W_R_conf", "gamma_ICU_W_R")
})

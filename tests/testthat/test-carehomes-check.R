context("carehomes (check)")

test_that("N_tot, N_tot2 and N_tot3 stay constant", {
  p <- carehomes_parameters(0, "uk")
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

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
  p <- carehomes_parameters(0, "england", beta_value = 0)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(19, 1, 101)))
})


test_that("everyone is infected when beta is large", {
  p <- carehomes_parameters(0, "england", beta_value = 1e9)
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(
    dust::dust_iterate(mod, seq(0, 400, by = 4))))
  expect_true(all(y$S[, 1, -1] == 0))
})


test_that("No one is infected if I and E are 0 at t = 0", {
  p <- carehomes_parameters(0, "england")
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y <- carehomes_initial(info, 1, p)$state
  y[info$index$I_asympt] <- 0
  y[info$index$R_pre] <- 0
  y[info$index$PCR_pos] <- 0

  mod$set_state(y)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(19, 1, 101)))
})


test_that("No one is hospitalised, no-one dies if p_sympt_ILI is 0", {
  p <- carehomes_parameters(0, "england")
  p$p_sympt_ILI[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$E > 0L))
  expect_true(all(y$I_ILI == 0))
  expect_true(all(y$I_hosp_R_unconf == 0))
  expect_true(all(y$I_hosp_R_conf == 0))
  expect_true(all(y$I_hosp_D_unconf == 0))
  expect_true(all(y$I_hosp_D_conf == 0))
  expect_true(all(y$I_ICU_R_unconf == 0))
  expect_true(all(y$I_ICU_R_conf == 0))
  expect_true(all(y$I_ICU_D_unconf == 0))
  expect_true(all(y$I_ICU_D_conf == 0))
  expect_true(all(y$I_triage_R_unconf == 0))
  expect_true(all(y$I_triage_R_conf == 0))
  expect_true(all(y$I_triage_D_unconf == 0))
  expect_true(all(y$I_triage_D_conf == 0))
  expect_true(all(y$R_stepdown_unconf == 0))
  expect_true(all(y$R_stepdown_conf == 0))
  expect_true(all(y$D_hosp == 0))
  expect_true(all(y$I_comm_D == 0))
  expect_true(all(y$D_comm == 0))
})


test_that("No one is hospitalised, no-one dies if psi_hosp_ILI is 0", {
  p <- carehomes_parameters(0, "england")
  p$psi_hosp_ILI[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$E > 0L))
  expect_true(any(y$I_ILI > 0))
  expect_true(all(y$I_hosp_R_unconf == 0))
  expect_true(all(y$I_hosp_R_conf == 0))
  expect_true(all(y$I_hosp_D_unconf == 0))
  expect_true(all(y$I_hosp_D_conf == 0))
  expect_true(all(y$I_ICU_R_unconf == 0))
  expect_true(all(y$I_ICU_R_conf == 0))
  expect_true(all(y$I_ICU_D_unconf == 0))
  expect_true(all(y$I_ICU_D_conf == 0))
  expect_true(all(y$I_triage_R_unconf == 0))
  expect_true(all(y$I_triage_R_conf == 0))
  expect_true(all(y$I_triage_D_unconf == 0))
  expect_true(all(y$I_triage_D_conf == 0))
  expect_true(all(y$R_stepdown_unconf == 0))
  expect_true(all(y$R_stepdown_conf == 0))
  expect_true(all(y$D_hosp == 0))
  expect_true(all(y$I_comm_D == 0))
  expect_true(all(y$D_comm == 0))
})


test_that("No one is hospitalised, no-one recovers in edge case", {
  p <- carehomes_parameters(0, "england")
  p$I0_asympt[] <- 0
  p$p_sympt_ILI[] <- 1
  p$p_hosp_ILI_step <- 1
  p$psi_hosp_ILI[] <- 1
  p$p_death_comm_step <- 1
  p$psi_death_comm[] <- 1
  p$p_asympt[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  ## Move initial infectives to ILI
  y0 <- carehomes_initial(info, 1, p)$state
  y0[info$index$I_ILI] <- y0[info$index$I_asympt]
  y0[info$index$I_asympt] <- 0

  mod$set_state(y0)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$I_ILI > 0))
  expect_true(all(y$I_hosp_R_unconf == 0))
  expect_true(all(y$I_hosp_R_conf == 0))
  expect_true(all(y$I_hosp_D_unconf == 0))
  expect_true(all(y$I_hosp_D_conf == 0))
  expect_true(all(y$I_ICU_R_unconf == 0))
  expect_true(all(y$I_ICU_R_conf == 0))
  expect_true(all(y$I_ICU_D_unconf == 0))
  expect_true(all(y$I_ICU_D_conf == 0))
  expect_true(all(y$I_triage_R_unconf == 0))
  expect_true(all(y$I_triage_R_conf == 0))
  expect_true(all(y$I_triage_D_unconf == 0))
  expect_true(all(y$I_triage_D_conf == 0))
  expect_true(all(y$R_stepdown_unconf == 0))
  expect_true(all(y$R_stepdown_conf == 0))
  expect_true(all(y$R == 0))
  expect_true(all(y$D_hosp == 0))
})


test_that("No one is hospitalised, no-one recovers in edge case 2", {
  p <- carehomes_parameters(0, "england")
  p$p_sympt_ILI[] <- 1
  p$p_hosp_ILI_step <- 1
  p$psi_hosp_ILI[] <- 1
  p$p_death_comm_step <- 1
  p$psi_death_comm[] <- 1
  p$p_asympt[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()

  ## Move initial infectives to ILI
  y0 <- carehomes_initial(info, 1, p)$state
  y0[info$index$I_ILI] <- y0[info$index$I_asympt]
  y0[info$index$I_asympt] <- 0

  mod$set_state(y0)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$I_ILI > 0))
  expect_true(all(y$I_hosp_R_unconf == 0))
  expect_true(all(y$I_hosp_R_conf == 0))
  expect_true(all(y$I_hosp_D_unconf == 0))
  expect_true(all(y$I_hosp_D_conf == 0))
  expect_true(all(y$I_ICU_R_unconf == 0))
  expect_true(all(y$I_ICU_R_conf == 0))
  expect_true(all(y$I_ICU_D_unconf == 0))
  expect_true(all(y$I_ICU_D_conf == 0))
  expect_true(all(y$I_triage_R_unconf == 0))
  expect_true(all(y$I_triage_R_conf == 0))
  expect_true(all(y$I_triage_D_unconf == 0))
  expect_true(all(y$I_triage_D_conf == 0))
  expect_true(all(y$R_stepdown_unconf == 0))
  expect_true(all(y$R_stepdown_conf == 0))
  expect_true(all(y$R == 0))
  expect_true(all(y$D_hosp == 0))
})


test_that("No one dies in the community if psi_death_comm is 0", {
  p <- carehomes_parameters(0, "england")
  p$psi_death_comm[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$I_ILI > 0))
  expect_true(all(y$I_comm_D == 0))
  expect_true(all(y$D_comm == 0))
})


test_that("forcing hospital route results in correct path", {
  ## We're going to try a number of very similar tests here, so a
  ## helper function will help run a model with given probabilities
  ## and verify that some compartments have cases (nonzero) and others
  ## area all zeros.
  helper <- function(prob_ICU_hosp, prob_death_ICU, prob_death_hosp_D,
                     expect_cases, expect_zero) {
    p <- carehomes_parameters(0, "england")
    p$p_ICU_hosp_step <- ifelse(is.null(prob_ICU_hosp),
                                    p$p_ICU_hosp_step, 1)
    p$psi_ICU_hosp[] <- prob_ICU_hosp %||% p$psi_ICU_hosp[]
    p$p_death_hosp_D_step <- ifelse(is.null(prob_death_hosp_D),
                                 p$p_death_hosp_D_step, 1)
    p$psi_death_hosp_D[] <- prob_death_hosp_D %||% p$psi_death_hosp_D[]
    p$p_death_ICU_step <- ifelse(is.null(prob_death_ICU),
                                 p$p_death_ICU_step, 1)
    p$psi_death_ICU[] <- prob_death_ICU %||% p$psi_death_ICU[]

    mod <- carehomes$new(p, 0, 1)
    info <- mod$info()
    mod$set_state(carehomes_initial(info, 1, p)$state)
    mod$set_index(integer(0))
    y <- mod$transform_variables(
      drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

    ## Save some work by using the total of confirmed and unconfirmed
    y$I_hosp_R <- y$I_hosp_R_unconf + y$I_hosp_R_conf
    y$I_hosp_D <- y$I_hosp_D_unconf + y$I_hosp_D_conf
    y$I_triage_R <- y$I_triage_R_unconf + y$I_triage_R_conf
    y$I_triage_D <- y$I_triage_D_unconf + y$I_triage_D_conf
    y$I_ICU_R <- y$I_ICU_R_unconf + y$I_ICU_R_conf
    y$I_ICU_D <- y$I_ICU_D_unconf + y$I_ICU_D_conf
    y$R_stepdown <- y$R_stepdown_unconf + y$R_stepdown_conf

    for (i in expect_cases) {
      expect_true(any(y[[i]] > 0), label = sprintf("Expected cases in %s", i))
    }
    for (i in expect_zero) {
      expect_true(all(y[[i]] == 0), label = sprintf("Expected zeros in %s", i))
    }
  }

  ## p_ICU_hosp = 0, p_death_hosp_D = 0 no-one goes into ICU, no deaths
  helper(0, NULL, 0, "I_hosp_R",
         c("I_hosp_D", "I_ICU_R", "I_ICU_D", "I_triage_R", "I_triage_D",
           "R_stepdown", "D_hosp"))

  ## p_death_hosp = 1, p_ICU_hosp = 0 no-one goes into ICU, no
  ## recovery in hospital
  helper(0, NULL, 1, "I_hosp_D",
         c("I_hosp_R", "I_ICU_R", "I_ICU_D", "I_triage_R", "I_triage_D",
           "R_stepdown"))

  ## p_death_ICU = 1, p_ICU_hosp = 1 no-one goes in hosp_D / hosp_R,
  ## no recovery from ICU
  helper(1, 1, NULL, "I_ICU_D",
         c("I_hosp_R", "I_hosp_D", "I_ICU_R", "R_stepdown"))

  ## p_death_ICU = 0, p_ICU_hosp = 1 no-one goes in hosp_D / hosp_R,
  ## no deaths
  helper(1, 0, NULL, "I_ICU_R",
         c("I_hosp_R", "I_hosp_D", "I_ICU_D", "D_hosp"))
})


test_that("No one seroconverts if p_seroconversion is 0", {
  p <- carehomes_parameters(0, "england")
  p$p_seroconversion[] <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$R_neg > 0))
  expect_true(all(y$R_pos == 0))
})


test_that("No one does not seroconvert and no one seroreverts
          if p_seroconversion is 1 and gamma_R_pos is 0", {
  p <- carehomes_parameters(0, "england")
  p$p_seroconversion[] <- 1
  ## set gamma_R_pos to 0 so no-one seroreverts
  p$gamma_R_pos <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(all(y$R_neg == 0))
  expect_true(any(y$R_pos > 0))
})


test_that("R_pre parameters work as expected", {
  helper <- function(p_R_pre_1, gamma_R_pre_1, gamma_R_pre_2) {
    p <- carehomes_parameters(0, "uk")
    p$p_R_pre_1 <- p_R_pre_1
    p$gamma_R_pre_1 <- gamma_R_pre_1
    p$gamma_R_pre_2 <- gamma_R_pre_2

    mod <- carehomes$new(p, 0, 1)
    info <- mod$info()

    y0 <- carehomes_initial(info, 1, p)$state
    y0[info$index$R_pre] <- 0

    mod$set_state(y0)
    mod$set_index(integer(0))
    mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))
  }

  ## p_R_pre = 1, expect no cases in R_pre_2 stream
  y <- helper(1, 1, 0.5)
  expect_true(all(y$R_pre[, 2, ] == 0))

  ## p_R_pre = 0, expect no cases in R_pre_1 stream
  y <- helper(0, 1, 0.5)
  expect_true(all(y$R_pre[, 1, ] == 0))

  ## gamma_R_pre_1 = gamma_R_pre_2 = 0, expect no cases in R_pos
  y <- helper(0.5, 0, 0)
  expect_true(all(y$R_pos == 0))
  expect_true(all(y$R_neg == 0))

  ## gamma_R_pre_1 = Inf, gamma_R_pre_2 = 0, expect progression in one
  ## time-step to R_neg/R_pos just from R_pre_1
  y <- helper(0.5, Inf, 0)
  n <- length(y$time)
  expect_equal(diff(t(apply(y$R_pos, c(1, 3), sum) + y$R_neg)),
               t(y$R_pre[, 1, -n]))

  ## gamma_R_pre_1 = 0, gamma_R_pre_2 = Inf, expect progression in one
  ## time-step to R_neg/R_pos just from R_pre_2
  y <- helper(0.5, 0, Inf)
  n <- length(y$time)
  expect_equal(diff(t(apply(y$R_pos, c(1, 3), sum) + y$R_neg)),
               t(y$R_pre[, 2, -n]))

  ## gamma_R_pre_1 = Inf, gamma_R_pre_2 = Inf, expect progression in
  ## one time-step to R_neg/R_pos from both R_pre_1 and R_pre_2
  y <- helper(0.5, Inf, Inf)
  n <- length(y$time)
  expect_equal(diff(t(apply(y$R_pos, c(1, 3), sum) + y$R_neg)),
               t(apply(y$R_pre[, , -n], c(1, 3), sum)))
})


test_that("setting a gamma to Inf results immediate progression", {
  helper <- function(gamma_name, progression_name, compartment_name) {
    p <- carehomes_parameters(0, "england")
    p[[gamma_name]] <- Inf
    p[[progression_name]] <- max(p[[progression_name]], 2)

    mod <- carehomes$new(p, 0, 1)

    info <- mod$info()
    mod$set_state(carehomes_initial(info, 1, p)$state, 0)
    mod$set_index(integer(0))
    y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

    y$I_hosp_R <- y$I_hosp_R_unconf + y$I_hosp_R_conf
    y$I_hosp_D <- y$I_hosp_D_unconf + y$I_hosp_D_conf
    y$I_triage_R <- y$I_triage_R_unconf + y$I_triage_R_conf
    y$I_triage_D <- y$I_triage_D_unconf + y$I_triage_D_conf
    y$I_ICU_R <- y$I_ICU_R_unconf + y$I_ICU_R_conf
    y$I_ICU_D <- y$I_ICU_D_unconf + y$I_ICU_D_conf
    y$R_stepdown <- y$R_stepdown_unconf + y$R_stepdown_conf

    z <- y[[compartment_name]]

    expect_true(any(z > 0))

    i <- seq_len(length(y$time) - 1L)
    if (length(dim(z)) == 4) {
      expect_equal(z[, 2, , i + 1], z[, 1, , i])
    } else {
      expect_equal(z[, 2, i + 1], z[, 1, i])
    }
  }

  helper("gamma_E", "s_E", "E")
  helper("gamma_asympt", "s_asympt", "I_asympt")
  helper("gamma_mild", "s_mild", "I_mild")
  helper("gamma_ILI", "s_ILI", "I_ILI")
  helper("gamma_triage", "s_triage", "I_triage_R")
  helper("gamma_triage", "s_triage", "I_triage_D")
  helper("gamma_hosp_R", "s_hosp_R", "I_hosp_R")
  helper("gamma_hosp_D", "s_hosp_D", "I_hosp_D")
  helper("gamma_ICU_R", "s_ICU_R", "I_ICU_R")
  helper("gamma_ICU_D", "s_ICU_D", "I_ICU_D")
  helper("gamma_comm_D", "s_comm_D", "I_comm_D")
  helper("gamma_stepdown", "s_stepdown", "R_stepdown")
  helper("gamma_R_pos", "s_R_pos", "R_pos")
  helper("gamma_PCR_pre", "s_PCR_pre", "PCR_pre")
  helper("gamma_PCR_pos", "s_PCR_pos", "PCR_pos")
})


test_that("setting a gamma to 0 results in no progression", {
  helper <- function(gamma_name, progression_name, compartment_name) {
    p <- carehomes_parameters(0, "england")
    p[[gamma_name]] <- 0
    p[[progression_name]] <- max(p[[progression_name]], 2)

    mod <- carehomes$new(p, 0, 1)

    info <- mod$info()
    mod$set_state(carehomes_initial(info, 1, p)$state, 0)
    mod$set_index(integer(0))
    y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

    y$I_hosp_R <- y$I_hosp_R_unconf + y$I_hosp_R_conf
    y$I_hosp_D <- y$I_hosp_D_unconf + y$I_hosp_D_conf
    y$I_triage_R <- y$I_triage_R_unconf + y$I_triage_R_conf
    y$I_triage_D <- y$I_triage_D_unconf + y$I_triage_D_conf
    y$I_ICU_R <- y$I_ICU_R_unconf + y$I_ICU_R_conf
    y$I_ICU_D <- y$I_ICU_D_unconf + y$I_ICU_D_conf
    y$R_stepdown <- y$R_stepdown_unconf + y$R_stepdown_conf

    z <- y[[compartment_name]]

    expect_true(any(z > 0))

    if (length(dim(z)) == 4) {
      expect_true(all(z[, 2, , ] == 0))
    } else {
      expect_true(all(z[, 2, ] == 0))
    }
  }

  p <- carehomes_parameters(0, "england")
  helper("gamma_E", "s_E", "E")
  helper("gamma_asympt", "s_asympt", "I_asympt")
  helper("gamma_mild", "s_mild", "I_mild")
  helper("gamma_ILI", "s_ILI", "I_ILI")
  helper("gamma_triage", "s_triage", "I_triage_R")
  helper("gamma_triage", "s_triage", "I_triage_D")
  helper("gamma_hosp_R", "s_hosp_R", "I_hosp_R")
  helper("gamma_hosp_D", "s_hosp_D", "I_hosp_D")
  helper("gamma_ICU_R", "s_ICU_R", "I_ICU_R")
  helper("gamma_ICU_D", "s_ICU_D", "I_ICU_D")
  helper("gamma_comm_D", "s_comm_D", "I_comm_D")
  helper("gamma_stepdown", "s_stepdown", "R_stepdown")
  helper("gamma_R_pos", "s_R_pos", "R_pos")
  helper("gamma_PCR_pre", "s_PCR_pre", "PCR_pre")
  helper("gamma_PCR_pos", "s_PCR_pos", "PCR_pos")
})


test_that("No one is unconfirmed, if p_admit_conf = 1", {
  p <- carehomes_parameters(0, "england")
  p$p_admit_conf_step <- 1
  p$psi_admit_conf[] <- 1
  p$gamma_triage <- Inf
  p$gamma_hosp_R <- Inf
  p$gamma_hosp_D <- Inf
  p$gamma_ICU_R <- Inf
  p$gamma_ICU_D <- Inf
  p$gamma_stepdown <- Inf

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state, 0)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(all(y$I_hosp_R_unconf == 0))
  expect_true(any(y$I_hosp_R_conf > 0))
  expect_true(all(y$I_hosp_D_unconf == 0))
  expect_true(any(y$I_hosp_D_conf > 0))
  expect_true(all(y$I_triage_R_unconf == 0))
  expect_true(any(y$I_triage_R_conf > 0))
  expect_true(all(y$I_triage_D_unconf == 0))
  expect_true(any(y$I_triage_D_conf > 0))
  expect_true(all(y$I_ICU_R_unconf == 0))
  expect_true(any(y$I_ICU_R_conf > 0))
  expect_true(all(y$I_ICU_D_unconf == 0))
  expect_true(any(y$I_ICU_D_conf > 0))
  expect_true(all(y$R_stepdown_unconf == 0))
  expect_true(any(y$R_stepdown_conf > 0))

  admit_conf <- apply(y$I_hosp_R_conf[, 1, , ] +
                      y$I_hosp_D_conf[, 1, , ] +
                      y$I_triage_R_conf[, 1, , ] +
                      y$I_triage_D_conf[, 1, , ], 1, sum)

  expect_true(all(diff(y$cum_admit_conf) == admit_conf[-1]))
  expect_true(all(y$cum_new_conf == 0))
})


test_that("No one is confirmed, if p_admit_conf = 0 and gamma_test = 0", {
  p <- carehomes_parameters(0, "england")
  p$psi_admit_conf[] <- 0
  p$gamma_test <- 0

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(carehomes_initial(info, 1, p)$state, 0)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$I_hosp_R_unconf > 0))
  expect_true(all(y$I_hosp_R_conf == 0))
  expect_true(any(y$I_hosp_D_unconf > 0))
  expect_true(all(y$I_hosp_D_conf == 0))
  expect_true(any(y$I_triage_R_unconf > 0))
  expect_true(all(y$I_triage_R_conf == 0))
  expect_true(any(y$I_triage_D_unconf > 0))
  expect_true(all(y$I_triage_D_conf == 0))
  expect_true(any(y$I_ICU_R_unconf > 0))
  expect_true(all(y$I_ICU_R_conf == 0))
  expect_true(any(y$I_ICU_D_unconf > 0))
  expect_true(all(y$I_ICU_D_conf == 0))
  expect_true(any(y$R_stepdown_unconf > 0))
  expect_true(all(y$R_stepdown_conf == 0))
  expect_true(all(y$admit_new_conf == 0))
  expect_true(all(y$cum_new_conf == 0))
})


test_that("Instant confirmation if p_admit_conf = 0 and gamma_test = Inf", {
  p <- carehomes_parameters(0, "england")
  p$psi_admit_conf[] <- 0
  p$gamma_test <- Inf
  p$gamma_triage <- Inf
  p$gamma_hosp_R <- Inf
  p$gamma_hosp_D <- Inf
  p$gamma_ICU_R <- Inf
  p$gamma_ICU_D <- Inf
  p$gamma_stepdown <- Inf

  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state

  ## We want to set I_ICU_R_unconf[, 1, ] and I_ICU_D_unconf[, 1, ] to 50
  y0[info$index$I_ICU_R_unconf[1:19]] <- 50
  y0[info$index$I_ICU_D_unconf[1:19]] <- 50

  mod$set_state(y0, 0)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))
  n <- length(y$time)

  ## Check hosp_R
  expect_true(all(y$I_hosp_R_conf[, 1, , ] == 0))
  expect_equal(y$I_hosp_R_conf[, 2, , -1], y$I_hosp_R_unconf[, 1, , -n])
  expect_true(all(y$I_hosp_R_unconf[, 2, , ] == 0))

  ## Check hosp_D
  expect_true(all(y$I_hosp_D_conf[, 1, , ] == 0))
  expect_equal(y$I_hosp_D_conf[, 2, , -1], y$I_hosp_D_unconf[, 1, , -n])
  expect_true(all(y$I_hosp_D_unconf[, 2, , ] == 0))

  ## Check triage_R
  expect_true(all(y$I_triage_R_conf[, 1, , ] == 0))
  expect_equal(y$I_triage_R_conf[, 2, , -1], y$I_triage_R_unconf[, 1, , -n])
  expect_true(all(y$I_triage_R_unconf[, 2, , ] == 0))

  ## Check triage_D
  expect_true(all(y$I_triage_D_conf[, 1, , ] == 0))
  expect_equal(y$I_triage_D_conf[, 2, , -1], y$I_triage_D_unconf[, 1, , -n])
  expect_true(all(y$I_triage_D_unconf[, 2, , ] == 0))

  ## Check ICU_R
  expect_equal(y$I_ICU_R_conf[, 2, , 2], y$I_ICU_R_unconf[, 1, , 1])
  expect_equal(y$I_ICU_R_conf[, 1, , -1], y$I_triage_R_conf[, 2, , -n])
  expect_true(all(y$I_ICU_R_unconf[, 2, , ] == 0))

  ## Check ICU_D
  expect_equal(y$I_ICU_D_conf[, 2, , 2], y$I_ICU_D_unconf[, 1, , 1])
  expect_equal(y$I_ICU_D_conf[, 1, , -1], y$I_triage_D_conf[, 2, , -n])
  expect_true(all(y$I_ICU_D_unconf[, 2, , ] == 0))

  ## Check stepdown
  expect_equal(y$R_stepdown_conf[, 2, 2], y$R_stepdown_unconf[, 1, 1])
  I_ICU_R_conf <- drop(y$I_ICU_R_conf)
  expect_equal(y$R_stepdown_conf[, 1, -1], I_ICU_R_conf[, 2, -n])
  expect_true(all(y$R_stepdown_unconf[, 2, ] == 0))

  new_conf <- apply(y$I_hosp_R_conf[, 2, , ] +
                    y$I_hosp_D_conf[, 2, , ] +
                    y$I_triage_R_conf[, 2, , ] +
                    y$I_triage_D_conf[, 2, , ], 2, sum)
  new_conf[2] <- new_conf[2] +
    sum(y$I_ICU_R_conf[, 2, , 2] +
        y$I_ICU_D_conf[, 2, , 2] +
        y$R_stepdown_conf[, 2, 2])
  expect_true(all(diff(y$cum_new_conf) == new_conf[-1]))

  expect_true(all(y$cum_admit_conf == 0))
})

test_that("tots all summed correctly ", {
  p <- carehomes_parameters(0, "england")
  mod <- carehomes$new(p, 0, 1)
  info <- mod$info()
  y0 <- carehomes_initial(info, 1, p)$state
  mod$set_state(carehomes_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))
  expect_true(all(y$general_tot == apply(y$I_triage_R_conf, 4, sum) +
                    apply(y$I_triage_D_conf, 4, sum) +
                    apply(y$I_hosp_R_conf, 4, sum) +
                    apply(y$I_hosp_D_conf, 4, sum) +
                    apply(y$R_stepdown_conf, 3, sum)))
  expect_true(all(y$I_ICU_tot == apply(y$I_ICU_R_conf, 4, sum) +
                    apply(y$I_ICU_D_conf, 4, sum)))
  expect_true(all(y$hosp_tot == y$I_ICU_tot + y$general_tot))
  expect_true(all(y$D_hosp_tot == apply(y$D_hosp, 2, sum)))
  expect_true(all(y$D_comm_tot == apply(y$D_comm, 2, sum)))
  expect_true(all(y$D_tot == y$D_hosp_tot + y$D_comm_tot))

  # check the positivity sums
  expect_true(all(y$sero_pos == apply(y$R_pos[4:13, , ], 3, sum)))
  expect_true(all(y$react_pos == apply(y$PCR_pos[2:18, , ], 3, sum)))
})

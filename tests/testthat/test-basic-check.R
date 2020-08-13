context("basic (check)")

test_that("N_tot stays constant", {
  p <- basic_parameters(0, "england")
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  n_tot <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$N_tot)
  expect_true(all(n_tot == sum(p$population)))
})


test_that("there are no infections when beta is 0", {
  p <- basic_parameters(0, "england", beta_value = 0)
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(17, 1, 101)))
})


test_that("everyone is infected when beta is very high", {
  p <- basic_parameters(0, "england", beta_value = 1e100)
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)
  expect_true(all(s[, , -1] == 0))
})


test_that("No one is infected if I and E are 0 at t = 0", {
  p <- basic_parameters(0, "england")
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  y <- basic_initial(info, 1, p)$state
  y[info$index$I_asympt] <- 0
  mod$set_state(y)
  mod$set_index(integer(0))
  s <- dust::dust_iterate(mod, seq(0, 400, by = 4), info$index$S)

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(17, 1, 101)))
})


test_that("No one is hospitalised if p_sympt_ILI is 0", {
  p <- basic_parameters(0, "england")
  p$p_sympt_ILI[] <- 0
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$E > 0))
  expect_true(all(y$I_ILI == 0))
  expect_true(all(y$I_hosp == 0))
  expect_true(all(y$I_ICU == 0))
  expect_true(all(y$R_hosp == 0))
  expect_true(all(y$D == 0))
})


test_that("No one goes to ICU and no deaths if p_recov_hosp is 1", {
  p <- basic_parameters(0, "england")
  p$p_recov_hosp[] <- 1
  p$p_death_hosp[] <- 0
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$I_hosp > 0))
  expect_true(all(y$I_ICU == 0))
  expect_true(all(y$R_hosp == 0))
  expect_true(all(y$D == 0))
})


test_that("p_death_hosp = 1, p_recov_hosp = 0: no icu, no recovery", {
  p <- basic_parameters(0, "england")
  p$p_recov_hosp[] <- 0
  p$p_death_hosp[] <- 1
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$I_hosp > 0))
  expect_true(all(y$I_ICU == 0))
  expect_true(all(y$R_hosp == 0))
  expect_true(any(y$D > 0))
})


test_that("if p_recov_ICU = 0, no-one recovers in hospital", {
  p <- basic_parameters(0, "england")
  p$p_recov_ICU[] <- 0
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(
    drop(dust::dust_iterate(mod, seq(0, 400, by = 4))))

  expect_true(any(y$I_hosp > 0))
  expect_true(all(y$R_hosp == 0))
})


test_that("if gamma_E is Inf, E cases must progress in 1 timestep", {
  p <- basic_parameters(0, "england")
  p$gamma_E <- Inf

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  ## NOTE: movement is in one *timestep* not one *day*, so can't use
  ## "by = 4" here to get daily output (and in similar tests below)
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$E)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$E > 0))
  expect_true(all(y$E[, 2, , j] == y$E[, 1, , i]))
})


test_that("if gamma_asympt is Inf, I_asympt must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_asympt <- Inf
  p$s_asympt <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$I_asympt)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_asympt > 0))
  expect_true(all(y$I_asympt[, 2, , j] == y$I_asympt[, 1, , i]))
})


test_that("if gamma_mild is Inf, I_mild cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_mild <- Inf
  p$s_mild <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$I_mild)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_mild > 0))
  expect_true(all(y$I_mild[, 2, , j] == y$I_mild[, 1, , i]))
})


test_that("if gamma_ILI is Inf, I_ILI cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_ILI <- Inf
  p$s_ILI <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$I_ILI)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_ILI > 0))
  expect_true(all(y$I_ILI[, 2, , j] == y$I_ILI[, 1, , i]))
})


test_that("if gamma_hosp is Inf, I_hosp cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_hosp <- Inf
  p$s_hosp <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$I_hosp)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_hosp > 0))
  expect_true(all(y$I_hosp[, 2, , j] == y$I_hosp[, 1, , i]))
})


test_that("if gamma_ICU is Inf, I_ICU cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_ICU <- Inf
  p$s_ICU <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$I_ICU)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_ICU > 0))
  expect_true(all(y$I_ICU[, 2, , j] == y$I_ICU[, 1, , i]))
})


test_that("if gamma_rec is Inf, R_hosp cases must progress in 1 time-step", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_rec <- Inf

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  i <- seq_len(dim(y$R_hosp)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$R_hosp > 0))
  expect_true(all(y$R_hosp[, 2, , j] == y$R_hosp[, 1, , i]))
})


test_that("if gamma_E is 0, E stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_E <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  ## NOTE: movement is in one *timestep* not one *day*, so can't use
  ## "by = 4" here to get daily output (and in similar tests below)
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$E[, 1, , ] > 0))
  expect_true(all(y$E[, 2, , ] == 0))
})


test_that("if gamma_asympt is 0, I_asympt stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_asympt <- 0
  p$s_asympt <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$I_asympt[, 1, , ] > 0))
  expect_true(all(y$I_asympt[, 2, , ] == 0))
})


test_that("if gamma_mild is 0, I_mild stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_mild <- 0
  p$s_mild <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$I_mild[, 1, , ] > 0))
  expect_true(all(y$I_mild[, 2, , ] == 0))
})


test_that("if gamma_ILI is 0, I_ILI stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_ILI <- 0
  p$s_ILI <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$I_ILI[, 1, , ] > 0))
  expect_true(all(y$I_ILI[, 2, , ] == 0))
})


test_that("if gamma_hosp is 0, I_hosp stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_hosp <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$I_hosp[, 1, , ] > 0))
  expect_true(all(y$I_hosp[, 2, , ] == 0))
})


test_that("if gamma_ICU is 0, I_ICU stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_ICU <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$I_ICU[, 1, , ] > 0))
  expect_true(all(y$I_ICU[, 2, , ] == 0))
})


test_that("if gamma_ICU is 0, I_ICU stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_rec <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(integer(0))
  y <- mod$transform_variables(drop(dust::dust_iterate(mod, 0:400)))

  expect_true(any(y$R_hosp[, 1, , ] > 0))
  expect_true(all(y$R_hosp[, 2, , ] == 0))
})

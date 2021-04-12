context("basic (check)")

test_that("N_tot stays constant", {
  p <- basic_parameters(0, "england")
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(info$index$N_tot)
  n_tot <- mod$simulate(seq(0, 400, by = 4))
  expect_vector_equal(n_tot, sum(p$population))
})


test_that("there are no infections when beta is 0", {
  p <- basic_parameters(0, "england", beta_value = 0)
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(17, 1, 101)))
})


test_that("everyone is infected when beta is very high", {
  p <- basic_parameters(0, "england", beta_value = 1e100)
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))
  expect_vector_equal(s[, , -1], 0)
})


test_that("No one is infected if I and E are 0 at t = 0", {
  p <- basic_parameters(0, "england")
  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  y <- basic_initial(info, 1, p)$state
  y[info$index$I_A] <- 0
  mod$set_state(y)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Susceptible population is never drawn down:
  expect_equal(s, array(s[, , 1], c(17, 1, 101)))
})


test_that("No one is hospitalised if p_C is 0", {
  p <- basic_parameters(0, "england")
  p$p_C[] <- 0
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$E > 0))
  expect_vector_equal(y$I_C, 0)
  expect_vector_equal(y$I_hosp, 0)
  expect_vector_equal(y$I_ICU, 0)
  expect_vector_equal(y$R_hosp, 0)
  expect_vector_equal(y$D, 0)
})


test_that("No one goes to ICU and no deaths if p_recov_hosp is 1", {
  p <- basic_parameters(0, "england")
  p$p_recov_hosp[] <- 1
  p$p_death_hosp[] <- 0
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    mod$simulate(seq(0, 400, by = 4)))

  expect_true(any(y$I_hosp > 0))
  expect_vector_equal(y$I_ICU, 0)
  expect_vector_equal(y$R_hosp, 0)
  expect_vector_equal(y$D, 0)
})


test_that("p_death_hosp = 1, p_recov_hosp = 0: no icu, no recovery", {
  p <- basic_parameters(0, "england")
  p$p_recov_hosp[] <- 0
  p$p_death_hosp[] <- 1
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_hosp > 0))
  expect_vector_equal(y$I_ICU, 0)
  expect_vector_equal(y$R_hosp, 0)
  expect_true(any(y$D > 0))
})


test_that("if p_recov_ICU = 0, no-one recovers in hospital", {
  p <- basic_parameters(0, "england")
  p$p_recov_ICU[] <- 0
  mod <- basic$new(p, 0, 1)

  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(
    drop(mod$simulate(seq(0, 400, by = 4))))

  expect_true(any(y$I_hosp > 0))
  expect_vector_equal(y$R_hosp, 0)
})


test_that("if gamma_E is Inf, E cases must progress in 1 timestep", {
  p <- basic_parameters(0, "england")
  p$gamma_E <- Inf

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  ## NOTE: movement is in one *timestep* not one *day*, so can't use
  ## "by = 4" here to get daily output (and in similar tests below)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  i <- seq_len(dim(y$E)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$E > 0))
  expect_vector_equal(y$E[, 2, , j], y$E[, 1, , i])
})


test_that("if gamma_A is Inf, I_A must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_A <- Inf
  p$k_A <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <-  mod$transform_variables(drop(mod$simulate(0:400)))

  i <- seq_len(dim(y$I_A)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_A > 0))
  expect_vector_equal(y$I_A[, 2, , j], y$I_A[, 1, , i])
})


test_that("if gamma_C is Inf, I_C cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_C <- Inf
  p$k_C <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  i <- seq_len(dim(y$I_C)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_C > 0))
  expect_vector_equal(y$I_C[, 2, , j], y$I_C[, 1, , i])
})


test_that("if gamma_hosp is Inf, I_hosp cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_hosp <- Inf
  p$k_hosp <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  i <- seq_len(dim(y$I_hosp)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_hosp > 0))
  expect_vector_equal(y$I_hosp[, 2, , j], y$I_hosp[, 1, , i])
})


test_that("if gamma_ICU is Inf, I_ICU cases must progress in 1 timestep", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_ICU <- Inf
  p$k_ICU <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  i <- seq_len(dim(y$I_ICU)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$I_ICU > 0))
  expect_vector_equal(y$I_ICU[, 2, , j], y$I_ICU[, 1, , i])
})


test_that("if gamma_rec is Inf, R_hosp cases must progress in 1 time-step", {
  ## This checks that progression groups work for these parameters,
  ## even though they are no longer the default
  p <- basic_parameters(0, "england")
  p$gamma_rec <- Inf

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  i <- seq_len(dim(y$R_hosp)[[4]] - 1)
  j <- i + 1L
  expect_true(any(y$R_hosp > 0))
  expect_vector_equal(y$R_hosp[, 2, , j], y$R_hosp[, 1, , i])
})


test_that("if gamma_E is 0, E stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_E <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  ## NOTE: movement is in one *timestep* not one *day*, so can't use
  ## "by = 4" here to get daily output (and in similar tests below)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$E[, 1, , ] > 0))
  expect_vector_equal(y$E[, 2, , ], 0)
})


test_that("if gamma_A is 0, I_A stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_A <- 0
  p$k_A <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$I_A[, 1, , ] > 0))
  expect_vector_equal(y$I_A[, 2, , ], 0)
})


test_that("if gamma_C is 0, I_C stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_C <- 0
  p$k_C <- 2

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$I_C[, 1, , ] > 0))
  expect_vector_equal(y$I_C[, 2, , ], 0)
})


test_that("if gamma_hosp is 0, I_hosp stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_hosp <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$I_hosp[, 1, , ] > 0))
  expect_vector_equal(y$I_hosp[, 2, , ], 0)
})


test_that("if gamma_ICU is 0, I_ICU stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_ICU <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$I_ICU[, 1, , ] > 0))
  expect_vector_equal(y$I_ICU[, 2, , ], 0)
})


test_that("if gamma_ICU is 0, I_ICU stay in progression stage 1", {
  p <- basic_parameters(0, "england")
  p$gamma_rec <- 0

  mod <- basic$new(p, 0, 1)
  info <- mod$info()
  mod$set_state(basic_initial(info, 1, p)$state)
  y <- mod$transform_variables(drop(mod$simulate(0:400)))

  expect_true(any(y$R_hosp[, 1, , ] > 0))
  expect_vector_equal(y$R_hosp[, 2, , ], 0)
})

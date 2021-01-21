context("carehomes (Rt)")

test_that("Can calculate Rt", {
  ## This does not really require that we have run the real model, as
  ## everything is deterministic. So here we'll generate some data
  ## from the model and use it so that we can compare against some
  ## exact numbers.
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  res <- carehomes_Rt(steps, y[, 1, ], p)
  expect_equal(res, d$outputs$rt_1)
  res_all <- carehomes_Rt_trajectories(steps, y, p)
  expect_equal(res_all, d$outputs$rt_all)

  ## Beta is returned, results correct length
  expect_identical(res$beta, rep(p$beta_step, length(steps)))
  expect_true(all(lengths(res) == length(steps)))

  ## Check the Rt calculation (from eff_Rt)
  expect_true(length(unique(res$Rt_all)) == 1)
  expect_true(length(unique(res$Rt_general)) == 1)

  ## Effective Rt lower than Rt
  expect_true(all(res$Rt_all >= res$eff_Rt_all))
  expect_true(all(res$Rt_general >= res$eff_Rt_general))

  ## General population effective Rt lower than total
  expect_true(all(res$eff_Rt_all >= res$eff_Rt_general))

  expect_equal(names(res), names(res_all))
  for (nm in names(res)) {
    expect_equal(res[[nm]], res_all[[nm]][, 1, drop = TRUE])
  }

  ## Date is returned
  expect_equal(res$date, res$step * p$dt)
})


test_that("validate inputs in rt calculation", {
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_Rt(steps, y[-1, 1, ], p),
    "Expected 'S' to have 19 rows = 19 groups x 1 vaccine classes",
    fixed = TRUE)
  expect_error(
    carehomes_Rt(steps, y[, 1, -1], p),
    "Expected 'S' to have 85 columns, following 'step'",
    fixed = TRUE)
})


test_that("validate inputs in rt trajectories calculation", {
  d <- reference_data_rt()

  p <- d$inputs$p
  steps <- d$inputs$steps
  y <- d$inputs$y

  expect_error(
    carehomes_Rt_trajectories(steps, y[, 1, ], p),
    "Expected a 3d array of 'S'",
    fixed = TRUE)
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , -1], p),
    "Expected 3rd dimension of 'S' to have length 85, following 'step'",
    fixed = TRUE)
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , ], list(p)),
    "Expected 2nd dimension of 'S' to have length 1, following 'pars'")
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , ], p, shared_parameters = FALSE),
    "If not using shared parameters, expected a unnamed list for 'pars'")
  expect_error(
    carehomes_Rt_trajectories(steps, y[, , ], list(p),
                              shared_parameters = TRUE),
    "If using shared parameters, expected a named list for 'pars'")
})


test_that("Can set initial time", {
  ## This test also checks that we can alter parameter inputs
  d <- reference_data_rt()

  steps <- d$inputs$steps
  y <- d$inputs$y

  step0 <- seq(steps[[1]], by = 1, length.out = ncol(y))

  p <- rep(list(d$inputs$p), ncol(y))
  for (i in seq_along(p)) {
    p[[i]]$initial_step <- step0[[i]]
  }

  res1 <- carehomes_Rt_trajectories(steps, y, p,
                                    initial_step_from_parameters = FALSE)
  expect_equal(res1$step, matrix(steps, length(steps), ncol(y)))

  res2 <- carehomes_Rt_trajectories(steps, y, p,
                                    initial_step_from_parameters = TRUE)
  expect_equal(res2$step[1, ], step0)
  expect_equal(res2$step[-1, ], res1$step[-1, ])
})


test_that("Can vary beta over time", {
  d <- reference_data_rt()

  steps <- d$inputs$steps
  y <- d$inputs$y

  p <- d$inputs$p
  dt <- p$dt
  initial_date <- p$initial_step * dt
  beta_date <- initial_date + c(0, 21, 62)
  beta_value <- p$beta_step * c(1, 0.5, 0.8)
  p$beta_step <- sircovid_parameters_beta(beta_date, beta_value, dt)
  p <- rep(list(p), ncol(y))
  p[[3]]$beta_step <- p[[3]]$beta_step / 2

  res <- carehomes_Rt_trajectories(steps, y, p,
                                   initial_step_from_parameters = FALSE)

  expect_equal(
    res$beta,
    cbind(sircovid_parameters_beta_expand(steps, p[[1]]$beta_step),
          sircovid_parameters_beta_expand(steps, p[[2]]$beta_step),
          sircovid_parameters_beta_expand(steps, p[[3]]$beta_step)))

  ## Check the Rt calculation (from eff_Rt) - compare the first test
  expect_true(length(unique(res$Rt_all)) > 1)
  expect_true(length(unique(res$Rt_general)) > 1)
  expect_true(all(res$Rt_all >= res$eff_Rt_all))
  expect_true(all(res$Rt_general >= res$eff_Rt_general))
})

test_that("Can compute Rt with larger timestep", {
  ### Currently this test is failing: 
  ### if we do not do the adjustment of parameters (using adjust_all_gammas)
  ### then the epidemic trajectories are similar but R values don't match
  ### if we do the adjustment of parameters (using adjust_all_gammas)
  ### the the R match but not the epidemic trajectories
  
  n <- c(1, 4, 100)
  pars <- lapply(n, function(n)
    carehomes_parameters(0, "london", steps_per_day = n))
  np <- 50 # use number large enough to get stable mean behaviour 
  
  ## adjust the gamma parameters to account for the discretisation
  pars <- lapply(pars, adjust_all_gammas)
  
  ## run model with different time steps
  mod <- lapply(pars, function(p) carehomes$new(p, 0, np, seed = 1L))
  end <- sircovid_date("2020-05-31")
  for (i in seq_along(n)) {
    info <- mod[[i]]$info()
    initial <- carehomes_initial(info, np, pars[[i]])
    mod[[i]]$set_state(initial$state, initial$step)
    mod[[i]]$set_index(carehomes_index(info)$run)
  }
  
  index <- mod[[1]]$info()$index$S # should be the same for both models
  f <- function(m, p) {
    steps <- seq(0, length.out = 100, by = 1 / p$dt)
    ## TODO: probably want to set a sensible index here:
    dust::dust_iterate(m, steps, index)
  }
  y <- Map(f, mod, pars)
  
  steps <-
    lapply(1:3, function(i) seq(0, length.out = 100, by = 1 / pars[[i]]$dt))
  
  rt_all_1 <- carehomes_Rt_trajectories(steps[[1]], y[[1]], pars[[1]])
  rt_all_2 <- carehomes_Rt_trajectories(steps[[2]], y[[2]], pars[[2]])
  rt_all_3 <- carehomes_Rt_trajectories(steps[[3]], y[[3]], pars[[3]])
  
  ## Below is some comparison to check that the outputs "agree" in
  ## the broadest sense as they're stochastic and *should not*
  ## disagree as we've made the integration much coarser.
  
  ## compare AR or rather n still susceptible
  mean_S_final_1 <- mean(apply(y[[1]][, , 100], 2, sum))
  mean_S_final_2 <- mean(apply(y[[2]][, , 100], 2, sum))
  mean_S_final_3 <- mean(apply(y[[3]][, , 100], 2, sum))
  # compute relative error using dt = 1/4 as reference
  rel_error_1 <- (mean_S_final_1 - mean_S_final_2) / mean_S_final_2 
  expect_true(rel_error_1 < 5e-2)
  rel_error_3 <- (mean_S_final_3 - mean_S_final_2) / mean_S_final_2 
  expect_true(rel_error_3 < 5e-2)
  
  ## compare mean reproduction numbers 
  mean_Rt_general_1 <- apply(rt_all_1$Rt_general, 1, mean) # with dt = 1
  mean_Rt_general_2 <- apply(rt_all_2$Rt_general, 1, mean) # with dt = 1/4
  mean_Rt_general_3 <- apply(rt_all_3$Rt_general, 1, mean) # with dt = 1/100
  # compute relative error using dt = 1/4 as reference as more precise 
  rel_error_1 <- (mean_Rt_general_1 - mean_Rt_general_2) / mean_Rt_general_2 
  expect_true(all(rel_error_1 < 1e-2))
  rel_error_3 <- (mean_Rt_general_3 - mean_Rt_general_2) / mean_Rt_general_2 
  expect_true(all(rel_error_3 < 1e-2))
  
  ## compare mean reproduction numbers 
  mean_Rt_all_1 <- apply(rt_all_1$Rt_all, 1, mean) # with dt = 1
  mean_Rt_all_2 <- apply(rt_all_2$Rt_all, 1, mean) # with dt = 1/4
  # compute relative error using dt = 1/4 as reference as more precise 
  rel_error <- (mean_Rt_all_1 - mean_Rt_all_2) / mean_Rt_all_2 
  expect_true(all(rel_error < 1e-2))
  
  ## compare mean reproduction numbers
  mean_eff_Rt_general_1 <- apply(rt_all_1$eff_Rt_general, 1, mean) # with dt = 1
  mean_eff_Rt_general_2 <- apply(rt_all_2$eff_Rt_general, 1, mean) # with dt = 1/4
  # compute relative error using dt = 1/4 as reference as more precise 
  rel_error <- (mean_eff_Rt_general_1 - mean_eff_Rt_general_2) / mean_eff_Rt_general_2 
  expect_true(all(rel_error < 1e-2))
  
  ## compare mean reproduction numbers
  mean_eff_Rt_all_1 <- apply(rt_all_1$eff_Rt_all, 1, mean) # with dt = 1
  mean_eff_Rt_all_2 <- apply(rt_all_2$eff_Rt_all, 1, mean) # with dt = 1/4
  # compute relative error using dt = 1/4 as reference as more precise 
  rel_error <- (mean_eff_Rt_all_1 - mean_eff_Rt_all_2) / mean_eff_Rt_all_2 
  expect_true(all(rel_error < 1e-2))
  
  # ## visual check
  # matplot(t(apply(rt_all_1$eff_Rt_all, 1, quantile, probs = c(0.025, 0.5, 0.975))),
  #         type = "l", lty = c(2, 1, 2), col = "black",
  #         xlab = "Time",
  #         ylab = "Effective Rt")
  # matplot(t(apply(rt_all_2$eff_Rt_all, 1, quantile, probs = c(0.025, 0.5, 0.975))),
  #         type = "l", lty = c(2, 1, 2), col = "blue", add = TRUE)

})

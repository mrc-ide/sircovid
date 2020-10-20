context("carehomes (vaccination)")

test_that("there are no infections if everyone is vaccinated with a vaccine
          preventing 100% of acquisition", {
  p <- carehomes_parameters(0, "england", rel_susceptibility = c(1, 0))
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
  expect_true(all(s[, 1, ] == 0))

  ## Noone changes compartment within the vaccinated individuals
  expect_true(all(s[, 2, ] == s[, 2, 1]))
})


test_that("Can calculate Rt with a vaccination class", {
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

  d <- reference_data_rt()

  skip("This test needs fixing!")
  rt_1 <- carehomes_Rt(steps, y[, 1, ], p)
  rt_all <- carehomes_Rt_trajectories(steps, y, p)

  expect_equal(rt_1, d$outputs$rt_1)
  expect_equal(rt_all, d$outputs$rt_all)
})

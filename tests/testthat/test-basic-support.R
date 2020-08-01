context("basic (support)")

test_that("basic progression parameters", {
  p <- basic_parameters_progression()
  expect_setequal(
    names(p),
    c("s_E", "s_asympt", "s_mild", "s_ILI", "s_hosp", "s_ICU", "s_rec",
      "gamma_E", "gamma_asympt", "gamma_mild", "gamma_ILI", "gamma_hosp",
      "gamma_ICU", "gamma_rec"))

  ## TODO: Lilith; you had said that there were some constraints
  ## evident in the fractional representation of these values - can
  ## you write tests here that reflect that?
})


test_that("basic_parameters returns a list of parameters", {
  date <- sircovid_date("2020-02-01")
  beta_date <- sircovid_date(c("2020-02-01", "2020-02-14", "2020-03-15"))
  beta_value <- c(3, 1, 2)

  p <- basic_parameters(date, "uk")
  expect_type(p, "list")
  expect_length(p$beta_step, 1)
  expect_identical(p$m, sircovid_transmission_matrix())

  progression <- basic_parameters_progression()
  expect_identical(p[names(progression)], progression)

  shared <- sircovid_parameters_shared(date, "uk", NULL, NULL)
  expect_identical(p[names(shared)], shared)
})


test_that("can tune the noise parameter", {
  p1 <- basic_parameters_observation()
  p2 <- basic_parameters_observation(1e4)
  expect_setequal(names(p1), names(p2))
  v <- setdiff(names(p1), "exp_noise")
  expect_mapequal(
})


test_that("basic_index identifies ICU and D_tot", {
  info <- list(time = 1, N_tot = 1, I_ICU_tot = 1, D_tot = 1, beta_out = 1,
               S = 17)
  expect_equal(basic_index(info), list(run = 3:4))
})


test_that("basic_index identifies ICU and D_tot in real model", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england")
  mod <- basic$new(p, 0, 10)
  info <- mod$info()
  index <- basic_index(info)
  expect_equal(index$run[[1]], which(names(info) == "I_ICU_tot"))
  expect_equal(index$run[[2]], which(names(info) == "D_tot"))
})


test_that("basic compare function returns NULL if no data present", {
  state <- c(0, 0) # ICU, D
  prev_state <- c(0, 0) # ICU, D
  observed <- list(itu = NA, deaths = NA)
  pars <- basic_parameters_observation(Inf)
  expect_null(basic_compare(state, prev_state, observed, pars))
})


test_that("observation function correctly combines likelihoods", {
  state <- rbind(10:15, 1:6) # ICU, D
  prev_state <- matrix(1, 2, 6)
  observed1 <- list(itu = 13, deaths = NA)
  observed2 <- list(itu = NA, deaths = 3)
  observed3 <- list(itu = 13, deaths = 3)

  pars <- basic_parameters_observation(Inf)
  ll1 <- basic_compare(state, prev_state, observed1, pars)
  ll2 <- basic_compare(state, prev_state, observed2, pars)
  ll3 <- basic_compare(state, prev_state, observed3, pars)

  expect_length(ll3, 6)

  ll_itu <- ll_nbinom(observed3$itu, pars$phi_ICU * state[1, ],
                      pars$k_ICU, Inf)
  ll_deaths <- ll_nbinom(observed3$deaths, pars$phi_death * (state[2, ] - 1),
                         pars$k_death, Inf)
  expect_equal(ll1, ll_itu)
  expect_equal(ll2, ll_deaths)
  expect_equal(ll3, ll1 + ll2)
})


test_that("can compute initial conditions", {
  p <- basic_parameters(sircovid_date("2020-02-07"), "england")
  mod <- basic$new(p, 0, 10)
  info <- mod$info()

  initial <- basic_initial(info, 10, p)
  expect_setequal(names(initial), c("state", "step"))
  expect_equal(initial$step, p$initial_step)

  ## TODO: this index faff simplifies away after odin.dust improvements
  len <- vnapply(info, prod)
  start <- cumsum(len) - len + 1L
  expect_equal(
    initial$state[[start[["N_tot"]]]],
    sum(p$population))

  i_S <- seq(start[["S"]], by = 1, length.out = info[["S"]])
  i_I <- seq(start[["I_asympt"]], by = 1, length.out = info[["I_asympt"]][[1]])
  expect_equal(
    initial$state[i_S] + initial$state[i_I],
    p$population)

  expect_equal(
    initial$state[i_I],
    append(rep(0, 16), 10, after = 3))
})

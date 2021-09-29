context("lancelot (IFR_t)")

test_that("Can calculate IFR_t", {
  ## This does not really require that we have run the real model, as
  ## everything is deterministic. So here we'll generate some data
  ## from the model and use it so that we can compare against some
  ## exact numbers.
  d <- reference_data_lancelot_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  S <- d$inputs$S
  I_weighted <- d$inputs$I_weighted

  res <- lancelot_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p)
  expect_equal(res, d$outputs$ifr_t_1)
  res_all <- lancelot_ifr_t_trajectories(steps, S, I_weighted, p)
  expect_equal(res_all, d$outputs$ifr_t_all)

  ## Results correct length
  expect_vector_equal(lengths(res), length(steps))

  expect_equal(names(res), names(res_all))
  for (nm in names(res)) {
    expect_equal(res[[nm]], res_all[[nm]][, 1, drop = TRUE])
  }

  ## no vaccination in the model so we expect with and without
  ## vaccination versions of IFR to be equal
  expect_equal(res$IFR_t_general, res$IFR_t_general_no_vacc)
  expect_equal(res$IHR_t_general, res$IHR_t_general_no_vacc)
  expect_equal(res$IFR_t_all, res$IFR_t_all_no_vacc)
  expect_equal(res$IHR_t_all, res$IHR_t_all_no_vacc)
  expect_equal(res_all$IFR_t_general, res_all$IFR_t_general_no_vacc)
  expect_equal(res_all$IHR_t_general, res_all$IHR_t_general_no_vacc)
  expect_equal(res_all$IFR_t_all, res_all$IFR_t_all_no_vacc)
  expect_equal(res_all$IHR_t_all, res_all$IHR_t_all_no_vacc)

  ## Date is returned
  expect_equal(res$date, res$step * p$dt)
})


test_that("validate inputs in ifr_t calculation", {
  d <- reference_data_lancelot_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  S <- d$inputs$S
  I_weighted <- d$inputs$I_weighted

  expect_error(
    lancelot_ifr_t(steps, S[-1, 1, ], I_weighted[, 1, ], p),
    "Expected 'S' to have 19 rows = 19 groups x 1 vacc classes",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t(steps, S[, 1, -1], I_weighted[, 1, ], p),
    "Expected 'S' to have 85 columns, following 'step'",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I_weighted[-1, 1, ], p),
    "Expected 'I_weighted' to have 19 rows",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I_weighted[, 1, -1], p),
    "Expected 'I_weighted' to have 85 columns, following 'step'",
    fixed = TRUE)
})


test_that("validate inputs in rt trajectories calculation", {
  d <- reference_data_lancelot_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  S <- d$inputs$S
  I_weighted <- d$inputs$I_weighted

  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, 1, ], I_weighted, p),
    "Expected a 3d array of 'S'",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I_weighted[, 1, ], p),
    "Expected a 3d array of 'I_weighted'",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, , -1], I_weighted, p),
    "Expected 3rd dim of 'S' to have length 85, given 'step'",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I_weighted[, , -1], p),
    "Expected 3rd dim of 'I_weighted' to have length 85, given 'step'",
    fixed = TRUE)
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, , ],
                                I_weighted[, 1, , drop = FALSE], list(p)),
    "Expected 2nd dim of 'S' to have length 1, given 'pars'")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, 1, , drop = FALSE],
                                I_weighted[, , ], list(p)),
    "Expected 2nd dim of 'I_weighted' to have length 1, given 'pars'")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, , ], I_weighted,
                                p, shared_parameters = FALSE),
    "If not using shared parameters, expected a unnamed list for 'pars'")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, , ], I_weighted, list(p),
                                shared_parameters = TRUE),
    "If using shared parameters, expected a named list for 'pars'")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S[, 1, , drop = FALSE], I_weighted, p,
                                shared_parameters = TRUE),
    "Expected 'S' and 'I_weighted' to have same length of 2nd dim")
})


test_that("Can set initial time", {
  ## This test also checks that we can alter parameter inputs
  d <- reference_data_lancelot_ifr_t()

  steps <- d$inputs$steps
  S <- d$inputs$S
  I_weighted <- d$inputs$I_weighted

  step0 <- seq(steps[[1]], by = 1, length.out = ncol(S))

  p <- rep(list(d$inputs$p), ncol(S))
  for (i in seq_along(p)) {
    p[[i]]$initial_step <- step0[[i]]
  }

  res1 <- lancelot_ifr_t_trajectories(steps, S, I_weighted, p,
                                      initial_step_from_parameters = FALSE)
  expect_equal(res1$step, matrix(steps, length(steps), ncol(S)))

  res2 <- lancelot_ifr_t_trajectories(steps, S, I_weighted, p,
                                      initial_step_from_parameters = TRUE)
  expect_equal(res2$step[1, ], step0)
  expect_equal(res2$step[-1, ], res1$step[-1, ])
})


test_that("can filter IFR_t to wanted types", {
  d <- reference_data_lancelot_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  S <- d$inputs$S
  I_weighted <- d$inputs$I_weighted

  expect_mapequal(
    lancelot_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p,
                   type = "IFR_t_general"),
    d$outputs$ifr_t_1[c("step", "date", "IFR_t_general")])
  expect_mapequal(
    lancelot_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p,
                   type = c("IFR_t_all", "IHR_t_all",
                            "IFR_t_general", "IHR_t_general")),
    d$outputs$ifr_t_1[c("step", "date", "IFR_t_all", "IHR_t_all",
                        "IFR_t_general", "IHR_t_general")])

  expect_mapequal(
    lancelot_ifr_t_trajectories(steps, S, I_weighted, p,
                                type = "IFR_t_general_no_vacc"),
    d$outputs$ifr_t_all[c("step", "date", "IFR_t_general_no_vacc")])
  expect_mapequal(
    lancelot_ifr_t_trajectories(steps, S, I_weighted, p,
                                type = c("IFR_t_all_no_vacc",
                                         "IHR_t_all_no_vacc",
                                         "IFR_t_general_no_vacc",
                                         "IHR_t_general_no_vacc")),
    d$outputs$ifr_t_all[c("step", "date", "IFR_t_all_no_vacc",
                          "IHR_t_all_no_vacc", "IFR_t_general_no_vacc",
                          "IHR_t_general_no_vacc")])
})


test_that("can't compute IFR_t for unknown types", {
  d <- reference_data_lancelot_ifr_t()

  p <- d$inputs$p
  steps <- d$inputs$steps
  S <- d$inputs$S
  I_weighted <- d$inputs$I_weighted

  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p,
                   type = "max_IFR_t_general"),
    "Unknown IFR/IHR type 'max_IFR_t_general', must match '")
  expect_error(
    lancelot_ifr_t_trajectories(steps, S, I_weighted, p,
                                type = "max_IFR_t_general"),
    "Unknown IFR/IHR type 'max_IFR_t_general', must match '")
  expect_error(
    lancelot_ifr_t(steps, S[, 1, ], I_weighted[, 1, ], p,
                   type = c("IFR_t_all", "ifr_t_general")),
    "Unknown IFR/IHR type 'ifr_t_general', must match '")
})


test_that("Parameters affect IFR/IHR/ALOS as expected", {

  ## Note that m_CHW and m_CHR have been changed from defaults to avoid
  ## having all care home residents infected
  p <- lancelot_parameters(sircovid_date("2020-02-07"), "england",
                           m_CHW = 3e-6, m_CHR = 3e-6)

  np <- 1L
  mod <- lancelot$new(p, 0, np, seed = 1L)

  initial <- lancelot_initial(mod$info(), 10, p)
  mod$set_state(initial$state, initial$step)
  index_S <- mod$info()$index$S
  index_I <- mod$info()$index$I_weighted


  end <- sircovid_date("2020-05-01") / p$dt
  steps <- seq(initial$step, end, by = 1 / p$dt)

  set.seed(1)
  y <- mod$simulate(steps)
  S <- y[index_S, , ]
  I <- y[index_I, , ]

  helper <- function(type, par_name, par_value_lower_ifr,
                     par_value_higher_ifr) {

    if (grepl("^p_", par_name)) {
      p[[par_name]][, ] <- par_value_lower_ifr
    } else {
      p[[par_name]] <- par_value_lower_ifr
    }
    ifr_lower <- lancelot_ifr_t(steps, S, I, p)

    if (grepl("^p_", par_name)) {
      p[[par_name]][, ] <- par_value_higher_ifr
    } else {
      p[[par_name]] <- par_value_higher_ifr
    }
    ifr_higher <- lancelot_ifr_t(steps, S, I, p)

    expect_vector_lt(ifr_lower[[type]], ifr_higher[[type]])
  }

  helper("IFR_t_all", "p_C_step", 0, 1)
  helper("IFR_t_all", "p_H_step", 0, 1)
  helper("IFR_t_all", "p_G_D_step", 0, 1)
  helper("IFR_t_all", "p_H_D_step", 0, 1)
  helper("IFR_t_all", "p_ICU_D_step", 0, 1)
  helper("IFR_t_all", "p_W_D_step", 0, 1)

  helper("IHR_t_all", "p_C_step", 0, 1)
  helper("IHR_t_all", "p_H_step", 0, 1)
  helper("IHR_t_all", "p_G_D_step", 1, 0)

  helper("ALOS", "gamma_ICU_pre_step", Inf, 1)
  helper("ALOS", "gamma_ICU_D_step", Inf, 1)
  helper("ALOS", "gamma_ICU_W_D_step", Inf, 1)
  helper("ALOS", "gamma_ICU_W_R_step", Inf, 1)
  helper("ALOS", "gamma_H_D_step", Inf, 1)
  helper("ALOS", "gamma_H_R_step", Inf, 1)
  helper("ALOS", "gamma_W_D_step", Inf, 1)
  helper("ALOS", "gamma_W_R_step", Inf, 1)

  helper("ALOS", "k_ICU_pre", 1, 2)
  helper("ALOS", "k_ICU_D", 1, 2)
  helper("ALOS", "k_ICU_W_D", 1, 2)
  helper("ALOS", "k_ICU_W_R", 1, 2)
  helper("ALOS", "k_H_D", 1, 2)
  helper("ALOS", "k_H_R", 1, 2)
  helper("ALOS", "k_W_D", 1, 2)
  helper("ALOS", "k_W_R", 1, 2)
})

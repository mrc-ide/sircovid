test_that("Single strain IFR excluding immunity calculated as expected", {

  helper <- function(p_C = NULL, p_H = NULL, p_G_D = NULL, p_ICU = NULL,
                     p_ICU_D = NULL, p_H_D = NULL, p_W_D = NULL) {

    ## include a 2nd vaccine class but it should have no impact
    p <- lancelot_parameters(0, "england",
                             rel_susceptibility = c(1, 0.5),
                             rel_p_sympt = c(1, 0.5),
                             rel_p_hosp_if_sympt = c(1, 0.5),
                             rel_p_death = c(1, 0.5),
                             waning_rate = 1 / 20)
    p$p_C_step[, ] <- p_C %||% p$p_C_step
    p$p_H_step[, ] <- p_H %||% p$p_H_step
    p$p_G_D_step[, ] <- p_G_D %||% p$p_G_D_step
    p$p_ICU_step[, ] <- p_ICU %||% p$p_ICU_step
    p$p_ICU_D_step[, ] <- p_ICU_D %||% p$p_ICU_D_step
    p$p_H_D_step[, ] <- p_H_D %||% p$p_H_D_step
    p$p_W_D_step[, ] <- p_W_D %||% p$p_W_D_step

    p
  }

  pars <- list()
  ## 1. Expect IFR = 0, IHR = 0, HFR = NA (no hospitalisations)
  pars[[1]] <- helper(p_C = 0)
  ## 2. Expect IFR = 0, IHR = 0, HFR = NA (no hospitalisations)
  pars[[2]] <- helper(p_C = 1, p_H = 0)
  ## 3. Expect IFR = 1, IHR = 0, HFR = NA (no hospitalisations)
  pars[[3]] <- helper(p_C = 1, p_H = 1, p_G_D = 1)
  ## 4. Expect IFR = 0, IHR = 1, HFR = 0
  pars[[4]] <- helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 0, p_H_D = 0)
  ## 5. Expect IFR = 1, IHR = 1, HFR = 1
  pars[[5]] <- helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 0, p_H_D = 1)
  ## 6. Expect IFR = 1, IHR = 1, HFR = 1
  pars[[6]] <-
    helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 1, p_ICU_D = 0, p_W_D = 1)
  ## 7. Expect IFR = 0, IHR = 1, HFR = 0
  pars[[7]] <-
    helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 1, p_ICU_D = 0, p_W_D = 0)
  ## 4. Expect IFR = 1, IHR = 1, HFR = 1
  pars[[8]] <-
    helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 1, p_ICU_D = 1)

  severity <- lancelot_ifr_excl_immunity(1, pars)

  expected_ifr <- c(0, 0, 1, 0, 1, 1, 0, 1)
  expected_ihr <- c(0, 0, 0, 1, 1, 1, 1, 1)
  expected_hfr <- c(NA, NA, NA, 0, 1, 1, 0, 1)

  expect_equal(severity$IFR[1, 1, ], expected_ifr)
  expect_equal(severity$IHR[1, 1, ], expected_ihr)
  expect_equal(severity$HFR[1, 1, ], expected_hfr)
})

test_that("Multistrain IFR excluding immunity calculated as expected", {

  helper <- function(strain_rel_p_sympt, strain_rel_p_hosp_if_sympt,
                     strain_rel_p_G_D, strain_rel_p_death, strain_rel_p_icu,
                     ifr_strain2, ihr_strain2, hfr_strain2) {

    ## include a 2nd vaccine class but it should have no impact
    p <- lancelot_parameters(0, "england",
                             strain_transmission = c(1, 1),
                             strain_rel_p_sympt = c(1, strain_rel_p_sympt),
                             strain_rel_p_hosp_if_sympt = c(1, strain_rel_p_hosp_if_sympt),
                             strain_rel_p_G_D = c(1, strain_rel_p_G_D),
                             strain_rel_p_death = c(1, strain_rel_p_death),
                             strain_rel_p_icu = c(1, strain_rel_p_icu),
                             rel_susceptibility = c(1, 0.5),
                             rel_p_sympt = c(1, 0.5),
                             rel_p_hosp_if_sympt = c(1, 0.5),
                             rel_p_death = c(1, 0.5),
                             waning_rate = 1 / 20)
    ## Fix values for strain 1 - resulting in IFR = 1, IHR = 0.5, HFR = 1
    p$p_C_step[, ] <- 1
    p$p_H_step[, ] <- 1
    p$p_G_D_step[, ] <- 0.5
    p$p_ICU_step[, ] <- 1
    p$p_ICU_D_step[, ] <- 1
    p$p_H_D_step[, ] <- 1
    p$p_W_D_step[, ] <- 1

    severity <- lancelot_ifr_excl_immunity(1, list(p))

    expect_equal(severity$IFR[1, , 1], c(1, ifr_strain2))
    expect_equal(severity$IHR[1, , 1], c(0.5, ihr_strain2))
    expect_equal(severity$HFR[1, , 1], c(1, hfr_strain2))

  }

  ## strain_rel_p_sympt = 0
  helper(0, NULL, NULL, NULL, NULL, 0, 0, NA)
  ## strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 0
  helper(1, 0, NULL, NULL, NULL, 0, 0, NA)
  ## strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1 &
  ## strain_rel_p_G_D = 2 (p_G_D for strain 2 is 1)
  helper(1, 1, 2, NULL, NULL, 1, 0, NA)
  ## strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1 &
  ## strain_rel_p_G_D = 0 & strain_rel_p_death = 0 & strain_rel_p_icu = 0
  helper(1, 1, 0, 0, 0, 0, 1, 0)
  ## strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1 &
  ## strain_rel_p_G_D = 0 & strain_rel_p_death = 1 & strain_rel_p_icu = 0
  helper(1, 1, 0, 1, 0, 1, 1, 1)
  ## strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1 &
  ## strain_rel_p_G_D = 0 & strain_rel_p_death = 0 & strain_rel_p_icu = 1
  helper(1, 1, 0, 0, 1, 0, 1, 0)
  ## strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1 &
  ## strain_rel_p_G_D = 0 & strain_rel_p_death = 1 & strain_rel_p_icu = 1
  helper(1, 1, 0, 1, 1, 1, 1, 1)

  ## Check probabilities are capped as expected
  helper(50, 50, 0, 50, 50, 1, 1, 1)
})

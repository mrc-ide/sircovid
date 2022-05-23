test_that("Single strain IFR excluding immunity calculated as expected", {

  helper <- function(p_C = NULL, p_H = NULL, p_G_D = NULL, p_ICU = NULL,
                     p_ICU_D = NULL, p_H_D = NULL, p_W_D = NULL) {

    # include a 2nd vaccine class but it should have no impact
    p <- lancelot_parameters(0, "england",
                             carehome_beds = 0,
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
  # 1. Expect IFR = 0, IHR = 0, HFR = NA (no hospitalisations)
  pars[[1]] <- helper(p_C = 0)
  # 2. Expect IFR = 0, IHR = 0, HFR = NA (no hospitalisations)
  pars[[2]] <- helper(p_C = 1, p_H = 0)
  # 3. Expect IFR = 1, IHR = 0, HFR = NA (no hospitalisations)
  pars[[3]] <- helper(p_C = 1, p_H = 1, p_G_D = 1)
  # 4. Expect IFR = 0, IHR = 1, HFR = 0
  pars[[4]] <- helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 0, p_H_D = 0)
  # 5. Expect IFR = 1, IHR = 1, HFR = 1
  pars[[5]] <- helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 0, p_H_D = 1)
  # 6. Expect IFR = 1, IHR = 1, HFR = 1
  pars[[6]] <-
    helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 1, p_ICU_D = 0, p_W_D = 1)
  # 7. Expect IFR = 0, IHR = 1, HFR = 0
  pars[[7]] <-
    helper(p_C = 1, p_H = 1, p_G_D = 0, p_ICU = 1, p_ICU_D = 0, p_W_D = 0)
  # 4. Expect IFR = 1, IHR = 1, HFR = 1
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

    # include a 2nd vaccine class but it should have no impact
    p <- lancelot_parameters(0, "england",
                             carehome_beds = 0,
                             strain_transmission = c(1, 1),
                             strain_rel_p_sympt = c(1, strain_rel_p_sympt),
                             strain_rel_p_hosp_if_sympt =
                               c(1, strain_rel_p_hosp_if_sympt),
                             strain_rel_p_G_D = c(1, strain_rel_p_G_D),
                             strain_rel_p_death = c(1, strain_rel_p_death),
                             strain_rel_p_icu = c(1, strain_rel_p_icu),
                             rel_susceptibility = c(1, 0.5),
                             rel_p_sympt = c(1, 0.5),
                             rel_p_hosp_if_sympt = c(1, 0.5),
                             rel_p_death = c(1, 0.5),
                             waning_rate = 1 / 20)
    # Fix values for strain 1 - resulting in IFR = 1, IHR = 0.5, HFR = 1
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

  # Set strain_rel_p_sympt = 0
  helper(0, NULL, NULL, NULL, NULL, 0, 0, NA)
  # Set strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 0
  helper(1, 0, NULL, NULL, NULL, 0, 0, NA)
  # Set strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1
  # & strain_rel_p_G_D = 2 (p_G_D for strain 2 is 1)
  helper(1, 1, 2, NULL, NULL, 1, 0, NA)
  # Set strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1
  # & strain_rel_p_G_D = 0 & strain_rel_p_death = 0 & strain_rel_p_icu = 0
  helper(1, 1, 0, 0, 0, 0, 1, 0)
  # Set strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1
  # & strain_rel_p_G_D = 0 & strain_rel_p_death = 1 & strain_rel_p_icu = 0
  helper(1, 1, 0, 1, 0, 1, 1, 1)
  # Set strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1
  # & strain_rel_p_G_D = 0 & strain_rel_p_death = 0 & strain_rel_p_icu = 1
  helper(1, 1, 0, 0, 1, 0, 1, 0)
  # Set strain_rel_p_sympt = 1 & strain_rel_p_hosp_if_sympt = 1
  # & strain_rel_p_G_D = 0 & strain_rel_p_death = 1 & strain_rel_p_icu = 1
  helper(1, 1, 0, 1, 1, 1, 1, 1)

  # Check probabilities are capped as expected
  helper(50, 50, 0, 50, 50, 1, 1, 1)
})

test_that("IFR excluding immunity outputted with correct dimensions", {

  # Dimensions of outputs should be n_steps x n_strains x n_pars

  p <- lancelot_parameters(1, "uk", carehome_beds = 0)
  severity <-  lancelot_ifr_excl_immunity(1, list(p))
  expect_equal(dim(severity$IFR), c(1, 1, 1))
  expect_equal(dim(severity$IHR), c(1, 1, 1))
  expect_equal(dim(severity$HFR), c(1, 1, 1))
  expect_equal(length(severity$step), 1)

  p <- lancelot_parameters(1, "uk", carehome_beds = 0)
  severity <-  lancelot_ifr_excl_immunity(1, list(p, p))
  expect_equal(dim(severity$IFR), c(1, 1, 2))
  expect_equal(dim(severity$IHR), c(1, 1, 2))
  expect_equal(dim(severity$HFR), c(1, 1, 2))
  expect_equal(length(severity$step), 1)

  p <- lancelot_parameters(1, "uk", carehome_beds = 0)
  severity <-  lancelot_ifr_excl_immunity(c(1, 3), list(p))
  expect_equal(dim(severity$IFR), c(2, 1, 1))
  expect_equal(dim(severity$IHR), c(2, 1, 1))
  expect_equal(dim(severity$HFR), c(2, 1, 1))
  expect_equal(length(severity$step), 2)

  p <- lancelot_parameters(1, "uk", carehome_beds = 0,
                           strain_transmission = c(1, 1))
  severity <-  lancelot_ifr_excl_immunity(1, list(p))
  expect_equal(dim(severity$IFR), c(1, 2, 1))
  expect_equal(dim(severity$IHR), c(1, 2, 1))
  expect_equal(dim(severity$HFR), c(1, 2, 1))
  expect_equal(length(severity$step), 1)

  p <- lancelot_parameters(1, "uk", carehome_beds = 0,
                           strain_transmission = c(1, 1))
  severity <-  lancelot_ifr_excl_immunity(c(1, 3, 5), list(p, p, p, p))
  expect_equal(dim(severity$IFR), c(3, 2, 4))
  expect_equal(dim(severity$IHR), c(3, 2, 4))
  expect_equal(dim(severity$HFR), c(3, 2, 4))
  expect_equal(length(severity$step), 3)

})


test_that("Can correctly calculate IFR excluding immunity at different steps", {

  p_severity <-
    lancelot_parameters_severity(0.25,
                                 p_H =
                                   list(value = c(0.4, 0.2),
                                        date = sircovid_date(c("2020-03-01",
                                                               "2020-04-01"))),
                                 p_ICU = list(value = 0),
                                 p_G_D = list(value = 0),
                                 p_H_D =
                                   list(value = c(0.5, 0.25),
                                        date = sircovid_date(c("2020-06-01",
                                                               "2020-07-01"))))

  p <- lancelot_parameters(1, "uk", carehome_beds = 0,
                           severity = p_severity)

  date <- sircovid_date(c("2020-02-01", "2020-05-01", "2020-08-01"))
  step <- date * 4

  severity <-  lancelot_ifr_excl_immunity(step, list(p))

  # IHR should halve between step[1] & step[2]
  expect_equal(severity$IHR[1, 1, 1], 2 * severity$IHR[2, 1, 1])
  expect_equal(severity$IHR[2, 1, 1], severity$IHR[3, 1, 1])
  # HFR should halve between step[2] & step[3]
  expect_equal(severity$HFR[1, 1, 1], severity$HFR[2, 1, 1])
  expect_equal(severity$HFR[2, 1, 1], 2 * severity$HFR[3, 1, 1])
  # IFR should halve between step[1] & step[2], and between step[2] & step[3]
  expect_equal(severity$IFR[1, 1, 1], 2 * severity$IFR[2, 1, 1])
  expect_equal(severity$IFR[2, 1, 1], 2 * severity$IFR[3, 1, 1])

})


test_that("If infections by age are weighted equally, IFR is a simple mean", {

  # Population even across age groups
  p <- lancelot_parameters(1, "uk", population = rep(500, 17),
                           carehome_beds = 0)
  # All mixing at the same rate
  p$m[1:17, 1:17] <- 0.1
  # probability of being symptomatic is the same across groups
  p$p_C_step[, ] <- p$p_C_step[1, 1]

  severity <- lancelot_ifr_excl_immunity(1, list(p))

  ihr <- p$p_C_step * p$p_H_step * (1 - p$p_G_D_step)
  hfr <- (1 - p$p_ICU_step) * p$p_H_D_step +
    p$p_ICU_step * (p$p_ICU_D_step + (1 - p$p_ICU_D_step) * p$p_W_D_step)
  ifr <- ihr * hfr + p$p_C_step * p$p_H_step * p$p_G_D_step

  expect_equal(severity$IHR[1, 1, 1], mean(ihr[1:17]))
  expect_equal(severity$IFR[1, 1, 1], mean(ifr[1:17]))
  # HFR gets weighted by IHR
  expect_equal(severity$HFR[1, 1, 1], weighted.mean(hfr[1:17], ihr[1:17]))

})

test_that("Validate IFR excluding immunity inputs", {
  p <- lancelot_parameters(1, "uk", carehome_beds = 0)
  expect_error(lancelot_ifr_excl_immunity(array(1, c(2, 2)), list(p)),
               "Expected 'step' to be a vector")

  p2 <- lancelot_parameters(1, "uk", carehome_beds = 0,
                            strain_transmission = c(1, 1))
  expect_error(lancelot_ifr_excl_immunity(1, list(p, p2)),
               "All parameter sets must have the same number of strains")

  p3 <- p
  p3$n_age_groups <- 5
  expect_error(lancelot_ifr_excl_immunity(1, list(p, p3)),
               "All parameter sets must have the same number of age groups")

  p4 <- p
  p$ICU_transmission <- 1
  expect_error(lancelot_ifr_excl_immunity(1, list(p)),
               "Cannot currently compute IFR if any of 'hosp_transmission'")
})


test_that("Severity by age is calculated parametrically", {

  helper <- function(p, i) {
    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    state <- lancelot_initial(info, 1, p)

    index_i <- array(info$index[[i]], info$dim[[i]])

    mod$update_state(state = state)
    mod$set_index(info$index[[i]])

    y <- mod$simulate(seq(0, 400, by = 4))

    expect_equal(length(y), prod(info$dim[[i]]) * 101)

    y <- array(y, c(info$dim[[i]], 101))
  }

  p <- lancelot_parameters(1, "uk", carehome_beds = 0)

  ihr <- p$p_C_step * p$p_H_step * (1 - p$p_G_D_step)
  hfr <- (1 - p$p_ICU_step) * p$p_H_D_step +
    p$p_ICU_step * (p$p_ICU_D_step + (1 - p$p_ICU_D_step) * p$p_W_D_step)
  ifr <- ihr * hfr + p$p_C_step * p$p_H_step * p$p_G_D_step

  # Severity is expected to be highest for the 80+
  y <- helper(p, "ifr_age")[c(1:17), 101]
  x <- ifr[which(!is.na(y))]
  expect_true(sum(y, na.rm = TRUE) == sum(x))

  y <- helper(p, "ihr_age")[c(1:17), 101]
  x <- ihr[which(!is.na(y))]
  expect_true(sum(y, na.rm = TRUE) == sum(x))

  y <- helper(p, "hfr_age")[c(1:17), 101]
  x <- hfr[which(!is.na(y))]
  expect_true(sum(y, na.rm = TRUE) == sum(x))

})


test_that("Multistrain severity has expected behaviour", {

  # Seeding parameters
  n_seeded_new_strain_inf <- 10
  start_date <- sircovid_date("2020-01-31")
  date_seeding <- start_date # seed both strains on same day

  # Helper function that runs model with p parameters
  helper <- function(p) {
    np <- 10
    mod <- lancelot$new(p, 0, np, seed = 1L)
    info <- mod$info()
    y0 <- lancelot_initial(info, 1, p)
    mod$update_state(state = y0)
    y <- mod$transform_variables(
      drop(mod$simulate(seq(0, 400, by = 4))))

    y
  }

  # Initial params - both strains have same severity
  p <- lancelot_parameters(start_date, "england",
                           strain_transmission = c(1, 1),
                           strain_seed_date = date_seeding,
                           strain_seed_size = n_seeded_new_strain_inf,
                           strain_seed_pattern = rep(1, 4))
  dim <- dim(p$rel_p_death)

  ## Severity for both strains is greater than 0
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] > 0, na.rm = TRUE))


  # Immunity vs strain 1 infection only
  p$rel_p_sympt <- array(matrix(c(0, 1, 1, 0),
                                nrow = dim[1], ncol = dim[2],
                                byrow = TRUE), dim = dim)
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] > 0, na.rm = TRUE))


  # Immunity vs strain 2 infection only
  p$rel_p_sympt <- array(matrix(c(1, 0, 0, 1),
                                nrow = dim[1], ncol = dim[2],
                                byrow = TRUE), dim = dim)
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] == 0, na.rm = TRUE))


  # Immunity vs both strains
  p$rel_p_sympt <- array(0, dim = dim)
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] == 0, na.rm = TRUE))


  # Immunity vs strain 1 hospitalisation only
  p$rel_p_sympt <- array(1, dim = dim)
  p$rel_p_hosp_if_sympt <- array(matrix(c(0, 1, 1, 0),
                                        nrow = dim[1], ncol = dim[2],
                                        byrow = TRUE), dim = dim)
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] > 0, na.rm = TRUE))


  # Immunity vs strain 2 hospitalisation only
  p$rel_p_hosp_if_sympt <- array(matrix(c(1, 0, 0, 1),
                                        nrow = dim[1], ncol = dim[2],
                                        byrow = TRUE), dim = dim)
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] == 0, na.rm = TRUE))


  # Immunity vs strain 1 death only
  p_deaths <- c("rel_p_death", "rel_p_G_D", "rel_p_H_D", "rel_p_W_D",
                "rel_p_ICU_D")
  p$rel_p_hosp_if_sympt <- array(1, dim = dim)
  for (i in p_deaths) {
    p[[i]] <- array(matrix(c(0, 1, 1, 0),
                           nrow = dim[1], ncol = dim[2],
                           byrow = TRUE), dim = dim)
  }
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] == 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] > 0, na.rm = TRUE))


  # Immunity vs strain 2 death only
  for (i in p_deaths) {
    p[[i]] <- array(matrix(c(1, 0, 0, 1),
                           nrow = dim[1], ncol = dim[2],
                           byrow = TRUE), dim = dim)
  }
  y <- helper(p)
  expect_true(all(y$ifr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ifr_strain[2, , -1] == 0, na.rm = TRUE))
  expect_true(all(y$ihr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$ihr_strain[2, , -1] > 0, na.rm = TRUE))
  expect_true(all(y$hfr_strain[1, , -1] > 0, na.rm = TRUE) &&
                all(y$hfr_strain[2, , -1] == 0, na.rm = TRUE))

})


test_that("Vaccination has expected behaviour against severity", {

  # Helper function that runs model with p parameters and
  # returns disag severity state i
  helper <- function(p, i) {
    mod <- lancelot$new(p, 0, 1, seed = 1L)
    info <- mod$info()
    state <- lancelot_initial(info, 1, p)

    index_i <- array(info$index[[i]], info$dim[[i]])

    mod$update_state(state = state)
    mod$set_index(info$index[[i]])

    y <- mod$simulate(seq(0, 400, by = 4))

    expect_equal(length(y), prod(info$dim[[i]]) * 101)

    y <- array(y, c(info$dim[[i]], 101))
  }

  # Initial params - no vaccine = no effect on severity
  p <- lancelot_parameters(0, "london",
                           rel_susceptibility = c(1, 1),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1))

  # Sanity check 1 - none changed class, all severity in class 1
  y <- helper(p, "ifr_disag")
  expect_true(all(y[, 2, ] == y[, 2, 1], na.rm = TRUE))
  expect_true(all(y[, 1, -1] > 0, na.rm = TRUE) &&
                all(y[, 2, -1] == 0, na.rm = TRUE))


  # Vaccine schedule
  region <- "london"
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)


  # Sanity check 2 - Everyone changed vacc_class, all severity in class 2
  p <- lancelot_parameters(0, region,
                           rel_susceptibility = c(1, 1),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)
  y <- helper(p, "ifr_disag")
  expect_true(all(y[, 1, -1] == 0, na.rm = TRUE) &&
                all(y[, 2, -1] > 0, na.rm = TRUE))


  ## Now for the proper tests - evaluation is done at the end of the model
  ## to allow for vaccination to take effect

  ## VE 99% vs susceptibility = lower IFR in class 2 than class 1
  vaccine_schedule <- test_vaccine_schedule(daily_doses = 5000,
                                            region = region,
                                            mean_days_between_doses = 1000,
                                            uptake = 1)

  p <- lancelot_parameters(0, region,
                           rel_susceptibility = c(1, 0.01),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)
  y <- helper(p, "ifr_disag")
  expect_true(sum(y[, 1, 101], na.rm = TRUE) > sum(y[, 2, 101], na.rm = TRUE))
  expect_true(all(y[, 1, 101] > 0, na.rm = TRUE))
  expect_true(all(y[, 2, 101] > 0, na.rm = TRUE))


  ## VE 100% vs death = no IFR/HFR in class 2, but IHR > 0 for both classes
  p <- lancelot_parameters(0, region,
                           rel_susceptibility = c(1, 1),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 0),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)
  y <- helper(p, "ifr_disag")
  expect_true(all(y[, 2, 101] == 0, na.rm = TRUE) &&
                all(y[, 1, 101] > 0, na.rm = TRUE))

  y <- helper(p, "ihr_disag")
  expect_true(all(y[, 2, 101] > 0, na.rm = TRUE) &&
                all(y[, 1, 101] > 0, na.rm = TRUE))

  y <- helper(p, "hfr_disag")
  expect_true(all(y[, 2, 101] == 0, na.rm = TRUE) &&
                all(y[, 1, 101] > 0, na.rm = TRUE))
})

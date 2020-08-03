context("carehomes (support)")

test_that("carehomes progression parameters", {
  p <- carehomes_parameters_progression()
  expect_setequal(
    names(p),
    c("s_E", "s_asympt", "s_mild", "s_ILI", "s_comm_D", "s_hosp_D",
      "s_hosp_R", "s_ICU_D", "s_ICU_R", "s_triage", "s_stepdown", "s_PCR_pos",
      "gamma_E", "gamma_asympt", "gamma_mild", "gamma_ILI", "gamma_comm_D",
      "gamma_hosp_D", "gamma_hosp_R", "gamma_ICU_D", "gamma_ICU_R",
      "gamma_triage", "gamma_stepdown", "gamma_R_pre_1", "gamma_R_pre_2",
      "gamma_test", "gamma_PCR_pos"))

  ## TODO: Lilith; you had said that there were some constraints
  ## evident in the fractional representation of these values - can
  ## you write tests here that reflect that?
})


test_that("carehomes_parameters returns a list of parameters", {
  date <- sircovid_date("2020-02-01")
  p <- carehomes_parameters(date, "uk")

  expect_type(p, "list")
  expect_equal(p$beta_step, 0.08)

  ## Transmission matrix is more complicated - see below
  expect_identical(p$m[1:17, 1:17], sircovid_transmission_matrix("uk"))
  expect_identical(
    p$m,
    carehomes_transmission_matrix(0.1, 4e-5, 5e-4, p$population))

  progression <- carehomes_parameters_progression()
  expect_identical(p[names(progression)], progression)

  shared <- sircovid_parameters_shared(date, "uk", NULL, NULL)
  ## NOTE: This is updated with CHR and CHW but may be renamed later;
  ## see comment in carehomes_parameters()
  shared$N_age <- 19L
  expect_identical(p[names(shared)], shared)

  severity <- carehomes_parameters_severity(NULL, p$population, 0.7)
  expect_identical(p[names(severity)], severity)

  expect_equal(
    p$observation,
    carehomes_parameters_observation(1e6, p$observation$N_tot_15_64))
  expect_equal(p$observation$N_tot_15_64, sum(p$N_tot[4:13]))

  extra <- setdiff(names(p),
                   c("m", "observation",
                     names(shared), names(progression), names(severity)))
  expect_setequal(
    extra,
    c("N_tot", "carehome_beds", "carehome_residents", "carehome_workers"))

  expect_equal(p$carehome_beds, sircovid_carehome_beds("uk"))
  expect_equal(p$carehome_residents, round(p$carehome_beds * 0.742))
  expect_equal(p$carehome_workers, p$carehome_residents)

  ## Can be slightly off due to rounding error
  expect_true(abs(sum(p$N_tot) - sum(p$population)) < 3)
  expect_length(p$N_tot, 19)
})


test_that("can compute severity for carehomes model", {
  population <- sircovid_population("uk")
  severity <- carehomes_parameters_severity(NULL, population, 0.7)
  expect_true(all(lengths(severity) == 19))
  expect_setequal(names(severity), names(sircovid_parameters_severity(NULL)))

  expect_true(
    all(severity$p_serocoversion == severity$p_serocoversion[[1]]))
  expect_equal(
    severity$p_death_comm, rep(c(0, 0.7), c(18, 1)))
  expect_equal(
    severity$p_admit_conf, rep(0.2, 19))
})


test_that("Can compute transmission matrix for carehomes model", {
  population <- sircovid_population("uk")
  m <- carehomes_transmission_matrix(0.1, 4e-5, 5e-4, population)
  expect_equal(rownames(m)[18:19], c("CHW", "CHR"))
  expect_equal(colnames(m)[18:19], c("CHW", "CHR"))
  expect_equal(dim(m), c(19, 19))
  expect_equal(m[1:17, 1:17], sircovid_transmission_matrix("uk"))
  expect_equal(unname(m[18:19, 18:19]),
               matrix(rep(c(4e-5, 5e-4), c(3, 1)), 2))
  expect_equal(m, t(m))

  ## TODO: we should get a check of the weighted mean and resident
  ## scaling here
})


test_that("can tune the noise parameter", {
  p1 <- carehomes_parameters_observation(1e6, 100)
  p2 <- carehomes_parameters_observation(1e4, 100)
  expect_setequal(names(p1), names(p2))
  v <- setdiff(names(p1), "exp_noise")
  expect_mapequal(p1[v], p2[v])
  expect_equal(p1$exp_noise, 1e6)
  expect_equal(p2$exp_noise, 1e4)
  expect_equal(p1$N_tot_15_64, 100)
})


test_that("carehomes_index identifies ICU and D_tot in real model", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 10)
  info <- mod$info()
  index <- carehomes_index(info)
  expect_equal(index$run[[1]], which(names(info$index) == "I_ICU_tot"))
  expect_equal(index$run[[2]], which(names(info$index) == "general_tot"))
  expect_equal(index$run[[3]], which(names(info$index) == "D_comm_tot"))
  expect_equal(index$run[[4]], which(names(info$index) == "D_hosp_tot"))
  expect_equal(index$run[[5]], which(names(info$index) == "D_tot"))
  expect_equal(index$run[[6]], which(names(info$index) == "cum_admit_conf"))
  expect_equal(index$run[[7]], which(names(info$index) == "cum_new_conf"))
  expect_equal(index$run[[8]], which(names(info$index) == "R_pre_15_64"))
  expect_equal(index$run[[9]], which(names(info$index) == "R_neg_15_64"))
  expect_equal(index$run[[10]], which(names(info$index) == "R_pos_15_64"))
})


test_that("carehome worker index is sensible", {
  expect_equal(carehomes_index_workers(), 6:13)
})


test_that("Can compute initial conditions", {
  p <- carehomes_parameters(sircovid_date("2020-02-07"), "england")
  mod <- carehomes$new(p, 0, 10)
  info <- mod$info()

  initial <- carehomes_initial(info, 10, p)
  expect_setequal(names(initial), c("state", "step"))
  expect_equal(initial$step, p$initial_step)

  initial_y <- mod$transform_variables(initial$state)

  expect_equal(initial_y$N_tot2, sum(p$population))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(initial_y$S + drop(initial_y$I_asympt), p$N_tot)
  expect_equal(drop(initial_y$I_asympt),
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$R_pre[, 1],
               append(rep(0, 18), 10, after = 3))
  expect_equal(initial_y$PCR_pos[, 1],
               append(rep(0, 18), 10, after = 3))

  ## 42 here, derived from;
  ## * 19 (S)
  ## * 19 (N_tot)
  ## * 4 values as N_tot2 + I_asympt[4] + R_pre[4] + PCR_pos[4]
  expect_equal(sum(initial$state != 0), 42)
})


test_that("get carehomes population can get the population", {
  expect_equal(sircovid_carehome_beds("uk"), 537521)
  expect_equal(sircovid_carehome_beds("UK"), 537521)
})


test_that("sircovid_carehome_beds throws sensible error on invalid input", {
  expect_error(sircovid_carehome_beds(NULL), "'region' must not be NULL")
  expect_error(
    sircovid_carehome_beds("oxfordshire"),
    "Carehome beds not found for 'oxfordshire': must be one of 'east_of_")
})


test_that("sircovid_carehome_beds caches data", {
  clear_cache()
  n <- sircovid_carehome_beds("uk")
  expect_s3_class(cache$carehomes, "data.frame")
  expect_equal(cache$carehomes$carehome_beds[cache$carehomes$region == "uk"],
               n)
})


## TODO: Ed - you had said that you had ideas for some more systematic
## testing here.  This will also get easier to do if we move to having
## a function generator given a data set.
test_that("carehomes_compare combines likelihood correctly", {
  state <- rbind(
    10:15, # ICU
    20:25, # general
    1:6,   # D_comm_tot
    3:8,   # D_hosp_tot
    4:9,   # D_tot
    50:55, # cum_admit_conf
    60:65, # cum_new_conf
    80:85, # R_pre_15_64
    30:35, # R_neg_15_64
    40:45) # R_pos_15_64
  prev_state <- array(1, dim(state))
  observed <- list(
    itu = 13,
    general = 23,
    deaths_hosp = 5,
    deaths_comm = 3,
    deaths = 8,
    admitted = 53,
    new = 63,
    npos_15_64 = 43,
    ntot_15_64 = 83)
  date <- sircovid_date("2020-01-01")
  pars <- carehomes_parameters(date, "uk", exp_noise = Inf)

  observed_keep <- function(nms) {
    observed[setdiff(names(observed), nms)] <- NA_real_
    observed
  }
  observed_drop <- function(nms) {
    observed[nms] <- NA_real_
    observed
  }

  ## This function is more complicated to test than the basic model
  ## because it's not a simple sum
  nms_sero <- c("npos_15_64", "ntot_15_64")
  parts <- c(as.list(setdiff(names(observed), nms_sero)), list(nms_sero))

  ll_parts <- lapply(parts, function(x)
    carehomes_compare(state, prev_state, observed_keep(x), pars))

  ## Extremely light testing, though this has already flushed out some
  ## issues
  expect_true(all(lengths(ll_parts) == 6))
  expect_equal(
    carehomes_compare(state, prev_state, observed, pars),
    rowSums(do.call(cbind, ll_parts)))
})


test_that("carehomes_population prevents negative populations", {
  population <- rep(100, 17)
  expect_error(
    carehomes_population(population, 1000, 0),
    "Not enough population to be care workers")
  expect_error(
    carehomes_population(population, 0, 1000),
    "Not enough population to meet care home occupancy")
})


test_that("carehomes_population preserves population", {
  population <- c(949, 989, 1030, 995, 993, 985, 965, 1082, 1042, 960,
                  980, 1004, 934, 1049, 971, 1020, 937)
  res <- carehomes_population(population, 100, 200)
  expect_equal(sum(res), sum(population))
  expect_equal(res[18:19], c(100, 200))
  expect_equal(res[1:5], population[1:5])
  expect_true(all(res[6:17] < population[6:17]))
})

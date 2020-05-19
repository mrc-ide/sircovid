context("parameters")

test_that("Parameters are generated as before", {
  time_steps_per_day <- 4

  pars_model <- generate_parameters(
    sircovid_model = basic_model(progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
                                 gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5)),
    transmission_model = "POLYMOD",
    severity_data_file = "extdata/severity_first.csv",
    beta = 0.125,
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/time_steps_per_day,
    use_polymod_pop = TRUE
  )

  # Remove newly generated parameters
  pars_model$beta_t <- NULL
  pars_model$beta_y <- NULL

  cmp <- readRDS("reference_pars.rds")

  # Remove parameters no longer generated
  pars_model$beta <- NULL
  pars_model$age_bin_starts <- NULL
  cmp$beta <- NULL
  
  expect_mapequal(pars_model, cmp)
})

test_that("Parameters generated as expected", {
  
  # Test adding transmission classes
  trans_profile = c(0.5, 0.2, 0.3)
  trans_increase = c(1, 2, 5)
  progression_groups = list(E = 3, asympt = 3, mild = 3, ILI = 3, hosp = 3, ICU = 3, rec = 3)
  
  test_params <- generate_parameters(sircovid_model = basic_model(progression_groups = progression_groups),
                                     trans_profile = trans_profile,
                                     trans_increase = trans_increase)
  expect_equal(test_params$trans_classes, length(trans_profile))
  expect_equal(dim(test_params$trans_profile), dim(test_params$trans_profile))
  
  expect_true(all(apply(test_params$trans_increase, 1, identical, trans_increase)))
  expect_true(all(apply(test_params$trans_profile, 1, identical, trans_profile)))
  
  # Check that partition dimensions are as expected
  expect_length(test_params$S0, test_params$N_age)
  expect_length(test_params$R0, test_params$N_age)
  expect_length(test_params$D0, test_params$N_age)
  expect_equal(dim(test_params$E0), c(test_params$N_age, progression_groups$E, length(trans_profile)))
  expect_equal(dim(test_params$I0_mild), c(test_params$N_age, progression_groups$mild, length(trans_profile)))
  expect_equal(dim(test_params$I0_asympt), c(test_params$N_age, progression_groups$asympt, length(trans_profile)))
  expect_equal(dim(test_params$I0_ILI), c(test_params$N_age, progression_groups$ILI, length(trans_profile)))
  expect_equal(dim(test_params$I0_hosp), c(test_params$N_age, progression_groups$hosp, length(trans_profile)))
  expect_equal(dim(test_params$I0_ICU), c(test_params$N_age, progression_groups$ICU, length(trans_profile)))
  expect_equal(dim(test_params$R0_hosp), c(test_params$N_age, progression_groups$rec, length(trans_profile)))
  
  # Test date interpolation, and starts from zero, accounts for leap years (2020)
  beta <- c(0.1, 0.2, 0.3)
  dt <- 0.25
  test_params <- generate_parameters(sircovid_model = basic_model(),
                                     beta=beta,
                                     beta_times=sircovid_date(c("2020-03-02", "2020-03-15", "2020-04-01")),
                                     dt = dt)
  expect_identical(test_params$beta_y, beta)
  expect_identical(test_params$beta_t, c(0, 13, 30)/dt)
  
  test_params <- generate_parameters(sircovid_model = basic_model(),
                                     beta=beta,
                                     beta_times=sircovid_date(c("2020-02-02", "2020-02-15", "2020-03-01")))
  expect_identical(test_params$beta_t, c(0, 13, 28)/dt)
})

test_that("Bad inputs", {
  expect_error(generate_parameters(trans_profile = c(0, 1)), 
               "Lengths of transmissibility class arguments mismatching")
  expect_error(generate_parameters(trans_profile = c(0.5, 1), trans_increase = c(1,2)), 
               "trans_profile proportions must sum to 1")
  expect_error(generate_parameters(beta_times = c(0, 1)), 
               "Length of beta mismatching with length of transition times")
  expect_error(generate_parameters(infection_seeding=list(values=c(1,10), bins=c('15 to 19'))), 
               "Each infection seeding value must correspond to one bin")
  expect_error(generate_parameters(infection_seeding=list(values=c(1,10), bins=c('15 to 19', '15 to 19'))), 
               "Seeding is into the same bin multiple times")
  expect_error(generate_parameters(transmission_model="OLYMOD"), 
               "Only POLYMOD transmission model implemented")
  expect_error(generate_parameters(sircovid_model = 
                                     basic_model(progression_groups=list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2))), 
               "progression_groups need to be defined for all partitions")
  expect_error(generate_parameters(sircovid_model = 
                                     basic_model(gammas=list(mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5))), 
               "gammas need to be defined for all partitions")
  expect_error(generate_parameters(beta = rep(0.1, 3),
                                   beta_times = sircovid_date(c("2020-02-02", "2021-03-01", "2020-04-01"))),
               "Supplied dates are not increasing")
})

test_that("Malformatted severity files are rejected", {
  expect_error(generate_parameters(severity_data_file="extdata/Final_COVID_severity.csv"), 
               "Could not find the following rows in the severity file:")
  expect_error(generate_parameters(severity_data_file="testdata/severity_test1.csv"), 
               "Not yet implemented decimal age bins")
  expect_message(generate_parameters(severity_data_file="testdata/severity_test2.csv"), 
               "Passed age bins intervals do not overlap correctly")
})

test_that("Time varying betas can be generated", {
  
  # Test time varying beta can be generated, and goes in and
  # out of generate_parameters() correctly
  beta = generate_beta(0.4, reduction_period = 3)
  dt <- 0.25
  params = generate_parameters(beta = beta$beta,
                               beta_times = beta$beta_times,
                               dt = dt)
  expect_equal(params$beta_t, c(0, 43, 44, 45)/dt)
  expect_equal(params$beta_y, c(0.4000, 0.4000, 0.2476, 0.0952))
  
  # Error checking
  expect_error(generate_beta(0.4, start_date = "2021-02-02", reduction_start = "2020-03-16"),
               "Start date must be earlier than intervention date")
  expect_error(generate_beta(-1),
               "beta_start must be non-negative")
  expect_error(generate_beta(0.4, beta_end = -1),
               "beta_end must be non-negative")
  expect_error(generate_beta(0.4, beta_reduction = -1),
               "beta cannot be reduced below 0")
  expect_message(generate_beta(0.4, reduction_period = 1000),
                 "Reduction period over 100 days - is this correct?")
    
})

test_that("read Bob's parameters", {
  ## just verify that the new system works without error so that we
  ## can continue with it.
  expect_error(
    generate_parameters(sircovid_model = hospital_model()),
    NA)
})

test_that("date conversion works", {
  first_data_date <- "2020-03-05"
  start_date <- "2020-03-01"
  
  offset <- start_date_to_offset(first_data_date, start_date)
  expect_equal(offset, 
               as.numeric(as.Date(first_data_date) - as.Date(start_date)))
  expect_equal(offset_to_start_date(first_data_date, offset),
               sircovid_date(start_date))
  expect_error(offset_to_start_date(first_data_date, start_date), 
               "Offset start date must be numeric")
  
})

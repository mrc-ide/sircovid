context("grid_search")

# Only tests that a grid search can be run
test_that("Small grid search works", {
  set.seed(1)

  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)

  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.2
  beta_step <- 0.05
  first_start_date <- "2020-01-21"
  last_start_date <- "2020-01-22"
  day_step <- 1
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  sircovid_model <- basic_model()
  model_params <- generate_parameters(
    sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = '2020-01-01',
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )
  
  scan_results <- scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    data = data,
    pars_obs = pars_obs,
    model_params = model_params,
    sircovid_model = sircovid_model,
  )

  expect_is(scan_results, "sircovid_scan")
  expect_true("inputs" %in% names(scan_results))
  expect_setequal(names(scan_results$inputs),
                  c("model", "model_params", "pars_obs", "data"))

  beta_grid <- seq(min_beta, max_beta, beta_step)
  date_grid <- seq(as.Date(first_start_date), as.Date(last_start_date), day_step)
  expect_equal(scan_results$x, beta_grid)
  expect_equal(scan_results$y, date_grid)
  expect_equal(dim(scan_results$renorm_mat_LL), dim(scan_results$mat_log_ll))
  expect_equal(dim(scan_results$renorm_mat_LL), c(length(beta_grid), length(date_grid)))
  expect_true(all(scan_results$renorm_mat_LL <= 1 & scan_results$renorm_mat_LL >= 0))

  # Plots run, but not checked
  plot(scan_results)
})

test_that("Warning is issued if grid does not explore low likelihood regions", {

  data <- readRDS("hospital_model_data.rds")
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  
  sircovid_model <- hospital_model()
  
  model_params <- generate_parameters(
    sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = '2020-01-01',
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )
  
  expect_warning(
    scan_beta_date(
    min_beta = 0.2,
    max_beta = 0.2,
    beta_step = 0.01,
    first_start_date = "2020-01-21",
    last_start_date = "2020-01-22",
    day_step = 1,
    data = data,
    sircovid_model = sircovid_model,
    pars_obs = pars_obs,
    model_params = model_params,
    scale_prior = 0.003541667,
    shape_prior = 36,
    n_particles = 2,
    ## 292 0s before a non-zero number when I run stuff.
    ## Hence I am setting tolerance to be very high.
    ## In practice, we migh need less.
    tolerance = 1e-10
    ),
    "Edges of the probability matrix are not close enough to 0."
  )
})


test_that("Small grid search works with new model", {
  set.seed(1)

  data <- readRDS("hospital_model_data.rds")

  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.2
  beta_step <- 0.05
  first_start_date <- "2020-01-21"
  last_start_date <- "2020-01-22"
  day_step <- 1
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  
  sircovid_model <- hospital_model()
  
  model_params <- generate_parameters(
    sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = '2020-01-01',
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )

  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    data = data,
    sircovid_model = sircovid_model,
    pars_obs = pars_obs,
    model_params = model_params,
    scale_prior = 0.003541667,
    shape_prior = 36
    )



  expect_is(scan_results, "sircovid_scan")
  expect_true("inputs" %in% names(scan_results))
  expect_setequal(names(scan_results$inputs),
                  c("model", "model_params", "pars_obs", "data"))

  beta_grid <- seq(min_beta, max_beta, beta_step)
  date_grid <- seq(as.Date(first_start_date), as.Date(last_start_date), day_step)
  expect_equal(scan_results$x, beta_grid)
  expect_equal(scan_results$y, date_grid)
  expect_equal(dim(scan_results$renorm_mat_LL), dim(scan_results$mat_log_ll))
  expect_equal(dim(scan_results$renorm_mat_LL), c(length(beta_grid), length(date_grid)))
  expect_true(all(scan_results$renorm_mat_LL <= 1 & scan_results$renorm_mat_LL >= 0))

  # Plots run, but not checked
  plot(scan_results, what = "likelihood")
  plot(scan_results, what = "probability")
})


test_that("Transmission is more likely", {
  set.seed(1)

  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)

  # Parameters for run
  min_beta <- 0
  max_beta <- 0.1
  beta_step <- 0.1
  first_start_date <- "2020-01-21"
  last_start_date <- "2020-01-21"
  day_step <- 1
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  
  sircovid_model <- basic_model()
  model_params <- generate_parameters(
    sircovid_model = sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = '2020-01-01',
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )

  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    pars_obs = pars_obs,
    model_params = model_params,
    data = data)

  # No transmission b = 0 much less likely than some transmission b = 0.1
  expect_lt(scan_results$renorm_mat_LL[[1]], scan_results$renorm_mat_LL[[2]])

})

test_that("Unreasonable start dates are less likely", {
  set.seed(1)

  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)

  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.1
  beta_step <- 0.1
  first_start_date <- "2020-01-01"
  last_start_date <- "2020-02-29"
  day_step <- 20
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  
  sircovid_model <- basic_model()
  model_params <- generate_parameters(
    sircovid_model = sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = '2020-01-01',
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )

  scan_results = scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    beta_step = beta_step,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    sircovid_model = sircovid_model,
    pars_obs = pars_obs,
    model_params = model_params,
    data = data)

  # Mid Jan start most likely
  expect_gt(scan_results$renorm_mat_LL[[2]], scan_results$renorm_mat_LL[[1]])
  expect_gt(scan_results$renorm_mat_LL[[2]], scan_results$renorm_mat_LL[[3]])

})

test_that("Bad parameters create errors", {
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)

  # Parameters for run
  min_beta <- 0.1
  max_beta <- 0.1
  beta_step <- 0.1
  first_start_date <- "2020-01-01"
  last_start_date <- "2020-02-29"
  day_step <- 20
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )


  sircovid_model <- basic_model(progression_groups = list(E = 2, asympt = 1, mild = 1, ILI = 1, hosp = 2, ICU = 2, rec = 2),
                                gammas = list(E = 1/2.5, asympt = 1/2.09, mild = 1/2.09, ILI = 1/4, hosp = 2/1, ICU = 2/5, rec = 2/5))
  model_params <- generate_parameters(
    sircovid_model = sircovid_model,
    transmission_model = "POLYMOD",
    hosp_transmission = 0,
    ICU_transmission = 0,
    beta = c(0.1, 0.1, 0.1),
    beta_times = c("2020-02-02", "2020-03-01", "2020-04-01"),
    trans_profile = 1,
    trans_increase = 1,
    dt = 0.25)

  expect_error(scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    model_params = model_params,
    pars_obs = pars_obs,
    beta_step = beta_step,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    data = data),
    "Set beta variation through generate_beta_func in sircovid_model, not model_params")

  model_params <- generate_parameters(
    sircovid_model = sircovid_model,
    transmission_model = "POLYMOD",
    hosp_transmission = 0,
    ICU_transmission = 0,
    beta = 0.1,
    beta_times = "2020-02-02",
    trans_profile = 1,
    trans_increase = 1,
    dt = 0.25)

  expect_error(scan_beta_date(
    min_beta = min_beta,
    max_beta = max_beta,
    model_params = model_params,
    pars_obs = pars_obs,
    beta_step = beta_step,
    first_start_date = first_start_date,
    last_start_date = last_start_date,
    day_step = day_step,
    data = data,
    scale_prior = 36),
    "If provided, both scale_prior and shape_prior must both be numeric")
})

context("Sampling and forecasts")

# Only tests that a grid search can be run
test_that("sample_grid_scan works", {
  set.seed(1)
  
  # grab the data
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  
  # filter and tidy it
  # Don't start the data until it is after thid date
  data <- data[data$date > as.Date("2020-02-29"),]
  
  # which dates should be set to NA
  na_dates <- as.Date(c("2020-03-01","2020-03-02","2020-03-03","2020-03-04",
                        "2020-03-05","2020-03-06","2020-03-07","2020-03-08",
                        "2020-03-15","2020-03-16","2020-03-19"))
  data$itu[data$date %in% na_dates] <- NA
  
  # bring the deaths back by 2 days as ECDC is out of sync
  data$deaths[head(seq_len(nrow(data)), -2)] <- tail(data$deaths, -2)
  data$deaths[tail(seq_len(nrow(data)), 2)] <- NA
  
  
  # Parameters for run
  min_beta <- 0.10
  max_beta <- 0.18
  beta_step <- 0.08
  first_start_date <- sircovid_date("2020-01-29")
  last_start_date <- sircovid_date("2020-02-14")
  day_step <- 12
  
  sircovid_model <- basic_model()

  
  model_params <- generate_parameters(
    sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = sircovid_date('2020-01-01'),
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
    pars_obs = list(phi_ICU = 0.95,
                    k_ICU = 2,
                    phi_death = 926 / 1019,
                    k_death = 2,
                    exp_noise = 1e6),
    model_params = model_params,
    sircovid_model = sircovid_model)

  n_sample_pairs <- 4 
  forecast_days <- 5
  res <- sample_grid_scan(scan_results = scan_results,
                          n_sample_pairs = n_sample_pairs, 
                          n_particles = 10,
                          forecast_days = forecast_days)

  model <- res$inputs$model$odin_model(user = res$inputs$model_params)
  # check length based on model and dates
  days_between <- length( min(res$param_grid$start_date) : tail(rownames(res$trajectories[,,1]),1))
  expect_equal(dim(res$trajectories), c(days_between, length(model$initial()), n_sample_pairs))
  
  # check the summary is as expected
  res_summary <- summary(res)
  # dates correct
  expect_equal(rownames(res_summary), 
               as.character(seq(from = sircovid_date(tail(data$date, 1)) + 1, 
                                to = sircovid_date(tail(data$date, 1)) + forecast_days, 
                                by = 1)))
  # quantiles increase
  expect_true(all(res_summary[1,c(-1,-2)] - head(res_summary[1,-1], -1) >= 0))
  
  ## Testing plotting is always a nightmare
  if (TRUE) {
    plot(res, what = "ICU")
    plot(res, what = "deaths")
  }
  
})


test_that("sample_grid_scan works with new model", {
  set.seed(1)
  
  # grab the data
  data <- readRDS("hospital_model_data.rds")
  
  
  # Parameters for run
  min_beta <- 0.10
  max_beta <- 0.18
  beta_step <- 0.08
  first_start_date <- sircovid_date("2020-01-29")
  last_start_date <- sircovid_date("2020-02-14")
  day_step <- 12
  
  sircovid_model <- hospital_model()
  
  model_params <- generate_parameters(
    sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = sircovid_date('2020-01-01'),
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
    pars_obs = list(
      phi_general = 0.95,
      k_general = 2,
      phi_ICU = 0.95,
      k_ICU = 2,
      phi_death = 926 / 1019,
      k_death = 2,
      exp_noise = 1e6
    ),
    model_params = model_params,
    sircovid_model = sircovid_model)
  
  n_sample_pairs <- 4 
  res <- sample_grid_scan(scan_results = scan_results,
                          n_sample_pairs = n_sample_pairs, 
                          n_particles = 10)
  
  model <- res$inputs$model$odin_model(user = res$inputs$model_params)
  # check length based on model and dates
  days_between <- length( min(res$param_grid$start_date) : tail(rownames(res$trajectories[,,1]),1))
  expect_equal(dim(res$trajectories), c(days_between, length(model$initial()), n_sample_pairs))
  
  ## Testing plotting is always a nightmare
  if (TRUE) {
    plot(res, what = "ICU")
    plot(res, what = "deaths")
    plot(res, what = "general")
  }
  
})


test_that("sample_pmcmc works with hospital model", {
  set.seed(1)
  
  # grab the data
  data <- readRDS("hospital_model_data.rds")
  
  sircovid_model <- hospital_model()
  model_params <- generate_parameters(
    sircovid_model = sircovid_model,
    severity_data_file = 'extdata/severity_2020_04_12.csv')
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  
  
  n_mcmc <- 10
  n_chains <- 2
  pars_to_sample <- data.frame(
    names=c('beta_start',
            'beta_end',
            'beta_pl', 
            'start_date',  
            'gamma_triage', 
            'gamma_hosp_R', 
            'gamma_hosp_D', 
            'gamma_ICU_R', 
            'gamma_ICU_D', 
            'gamma_stepdown'),
    init=c(0.14, 
           0.14*0.238,
           0.14*0.238,
           sircovid_date("2020-02-07"),
           0.5099579,
           0.1092046,
           0.2911154,
           0.3541429,
           0.2913861,
           0.452381),
    min=c(0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0),
    max=c(1,
          1,
          1,
          1e6,
          1,
          1,
          1,
          1,
          1,
          1),
    discrete=c(FALSE,
               FALSE,
               FALSE,
               TRUE,
               FALSE,
               FALSE,
               FALSE,
               FALSE,
               FALSE,
               FALSE),
    stringsAsFactors = FALSE)
  pars_lprior <- list('beta_start'     = function(pars) log(1e-10),
                      'beta_end'       = function(pars) 0,
                      'beta_pl'        = function(pars) 0,
                      'start_date'     = function(pars) 0,
                      'gamma_triage'   = function(pars) 0,
                      'gamma_hosp_R'   = function(pars) 0,
                      'gamma_hosp_D'   = function(pars) 0,
                      'gamma_ICU_R'    = function(pars) 0,
                      'gamma_ICU_D'    = function(pars) 0,
                      'gamma_stepdown' = function(pars) 0)
  proposal_kernel <- diag(length(pars_to_sample$names)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- pars_to_sample$names
  proposal_kernel['start_date', 'start_date'] <- 25

  mcmc_results <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    pars_to_sample = pars_to_sample,
    pars_lprior = pars_lprior,
    proposal_kernel = proposal_kernel,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    n_chains = n_chains
  )

  
  n_sample <- 2
  res <- sample_pmcmc(mcmc_results = mcmc_results,
                      burn_in = 1,
                      n_sample = n_sample, 
                      n_particles = 10,
                      forecast_days = 0)
  
  model <- res$inputs$model$odin_model(user = res$inputs$model_params)
  # check length based on model and dates
  days_between <- length( sircovid_date(min(res$param_grid$start_date)) : tail(rownames(res$trajectories[,,1]),1))
  expect_equal(dim(res$trajectories), c(days_between, length(model$initial()), n_sample))
  
  # check forecasting
  forecast_days <- 2
  res <- sample_pmcmc(mcmc_results = mcmc_results,
                      burn_in = 1,
                      n_sample = n_sample, 
                      n_particles = 10,
                      forecast_days = forecast_days)
  expected_total_days <- tail(sircovid_date(data$date), 1) + forecast_days - sircovid_date(min(res$param_grid$start_date)) + 1
  expect_equal(dim(res$trajectories)[1], expected_total_days)
  
  ## Testing plotting
  if (TRUE) {
    plot(res, what = "ICU")
    plot(res, what = "deaths")
    plot(res, what = "general")
  }
  
})


test_that("sample_pmcmc works with serology model", {

  
  # grab the data
  data <- readRDS("serology_model_data.rds")
  sircovid_model <- serology_model()
  model_params <- generate_parameters(
    sircovid_model = sircovid_model,
    transmission_model = "POLYMOD",
    beta = 0.1,
    beta_times = sircovid_date('2020-01-01'),
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )
  pars_obs <-  list(phi_general = 0.95,
                    k_general = 2,
                    phi_ICU = 0.95,
                    k_ICU = 2,
                    phi_death = 1.15,
                    k_death = 2,
                    phi_new = 0.95,
                    k_new = 2,
                    phi_admitted = 0.95,
                    k_admitted = 2,
                    exp_noise = 1e6)
  
  par_names <- c('beta_start',
                 'beta_end', 
                 'beta_pl',
                 'start_date',  
                 'gamma_triage', 
                 'gamma_hosp_R', 
                 'gamma_hosp_D', 
                 'gamma_ICU_R', 
                 'gamma_ICU_D', 
                 'gamma_stepdown')
  
  
  n_mcmc <- 10
  n_particles <- 10
  set.seed(1)
  proposal_kernel <- diag(10)*0.001^2
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- par_names
  proposal_kernel["start_date", "start_date"] <- 1
  
  
  mcmc_results <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    proposal_kernel = proposal_kernel, 
    n_particles = n_particles, 
    n_chains = 2
  )

  
  n_sample <- 2
  res <- sample_pmcmc(mcmc_results = mcmc_results,
                      burn_in = 1,
                      n_sample = n_sample, 
                      n_particles = 10,
                      forecast_days = 0)
  
  model <- res$inputs$model$odin_model(user = res$inputs$model_params)
  # check length based on model and dates
  days_between <- length( sircovid_date(min(res$param_grid$start_date)) : tail(rownames(res$trajectories[,,1]),1))
  expect_equal(dim(res$trajectories), c(days_between, length(model$initial()), n_sample))
  
  # check forecasting
  forecast_days <- 2
  res <- sample_pmcmc(mcmc_results = mcmc_results,
                      burn_in = 1,
                      n_sample = n_sample, 
                      n_particles = 10,
                      forecast_days = forecast_days)
  expected_total_days <- tail(sircovid_date(data$date), 1) + forecast_days - sircovid_date(min(res$param_grid$start_date)) + 1
  expect_equal(dim(res$trajectories)[1], expected_total_days)
  
  ## Testing summary
  summary(res, what = "deaths")
  summary(res, what = "icu")
  summary(res, what = "hosp")

  ## Testing plotting
  if (TRUE) {
    plot(res, what = "ICU")
    plot(res, what = "deaths")
    plot(res, what = "general")
    plot(res, what = "admitted")
    plot(res, what = "new")
  }
  
})

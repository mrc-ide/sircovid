context("pmcmc")

test_that("pmcmc runs with beta_pl", {
  
  data <- readRDS("hospital_model_data.rds")
  sircovid_model <- hospital_model()
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

  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'beta_pl', 'start_date'),
    init=c(0.14, 0.14*0.238, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0, 0),
    max=c(1, 1, 1, 1e6),
    discrete=c(FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'beta_end' = function(pars) 0,
                     'beta_pl' = function(pars) 0,
                     'start_date' = function(pars) 0)

  proposal_kernel <- diag(nrow(pars_to_sample)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- pars_to_sample$names
  proposal_kernel['start_date', 'start_date'] <- 25
  
  set.seed(2)
  X2 <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    pars_to_sample = pars_to_sample,
    proposal_kernel = proposal_kernel,
    pars_lprior = pars_lprior,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs
  )
  
  expect_is(X2, 'pmcmc')
  expect_setequal(names(X2), c('inputs', 'results', 'states', 'acceptance_rate', 'ess'))
  expect_equal(dim(X2$results), c(n_mcmc + 1L, nrow(pars_to_sample) + 3L))
  expect_equal(dim(X2$states), c(n_mcmc + 1L, 289))
  
  
})

test_that("pmcmc runs without error", {
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)

  cmp <- readRDS("reference_pmcmc.rds")

  n_mcmc <- 10
  set.seed(1)
  sircovid_model <- basic_model(gammas = list(E = 1/(2.5),
                                              asympt = 1/2.09,
                                              mild = 1/2.09,
                                              ILI = 1/4,
                                              hosp = 2,
                                              ICU = 2/5,
                                              rec = 2/5))
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
  
  pars_obs <- list(phi_ICU = 0.95,
                   k_ICU = 2,
                   phi_death = 926 / 1019,
                   k_death = 2,
                   exp_noise = 1e6)

  pars_to_sample <- data.frame(
    names=c('beta_start', 'start_date'),
    init=c(0.14, sircovid_date("2020-02-07")),
    min=c(0, 0),
    max=c(1, 1e6),
    discrete=c(FALSE, TRUE),
    stringsAsFactors = FALSE)
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'start_date' = function(pars) 0)
  
  X <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    pars_to_sample = pars_to_sample,
    proposal_kernel = matrix(c(0.001^2, 0,
                       0, 0.5^2),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(
                       c('beta_start', 'start_date'),
                       c('beta_start', 'start_date'))),
    pars_lprior = pars_lprior
  )
  
  expect_is(X, 'pmcmc')
  expect_setequal(names(X), c('inputs', 'results', 'states', 'acceptance_rate', 'ess'))
  expect_equal(dim(X$results), c(n_mcmc + 1L, 5))
  expect_equal(dim(X$states), c(n_mcmc + 1L, 238))
  expect_equivalent(X[-1], cmp[-1])
  
  # Plots run, but not checked
  plot(X)
  # run summary method
  summary(X)
  
  
  cmp <- readRDS('reference_pmcmc_three_par.rds')

  ## test three par version

  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'beta_end' = function(pars) 0,
                     'start_date' = function(pars) 0)
  
  
 set.seed(2)
 X2 <- pmcmc(
   data = data,
   n_mcmc = n_mcmc,
   pars_to_sample = pars_to_sample,
   proposal_kernel = matrix(c(0.001^2, 0, 0,
                              0, 0.001^2, 0,
                              0,       0, 0.5^2), 
                            nrow = 3, byrow = TRUE,
                            dimnames = list(
                              pars_to_sample$names,
                              pars_to_sample$names)),
   pars_obs = pars_obs, 
   pars_lprior = pars_lprior,
   model_params = model_params,
   sircovid_model = sircovid_model
 )

 expect_equal(dim(X2$results), c(n_mcmc + 1L, 6))
 expect_equivalent(X2[-1], cmp[-1])

 ## check that proposing jumps of size zero results in the initial parameter being retained
 Z <- pmcmc(
   data = data,
   n_mcmc = n_mcmc,
   pars_to_sample = pars_to_sample,
   proposal_kernel = matrix(rep(0, 9), 
                            nrow = 3, byrow = TRUE,
                            dimnames = list(
                              pars_to_sample$names,
                              pars_to_sample$names)),
   pars_obs = pars_obs, 
   pars_lprior = pars_lprior,
   model_params = model_params,
   sircovid_model = sircovid_model,
   output_proposals = TRUE
 )

 expect_equal(object = Z$results$beta_start,
              expected = rep(Z$inputs$pars$pars_init$beta_start, n_mcmc + 1))
 expect_equal(object = Z$results$start_date,
              expected = rep(Z$inputs$pars$pars_init$start_date, n_mcmc + 1))
 expect_true(!all(diff(Z$results$log_likelihood) == 0))
 expect_equal(dim(Z$proposals), c(n_mcmc + 1L, 7))

## check non-zero covariance ihas an impact on proposals
set.seed(2)
X3 <- pmcmc(
  data = data,
  n_mcmc = n_mcmc,
  pars_to_sample = pars_to_sample,
  proposal_kernel = matrix(c(0.001^2, 0, 0.001*0.5*0.6,
                             0, 0.001^2, 0,
                             0.001*0.5*0.6, 0, 0.5^2), 
                           nrow = 3, byrow = TRUE,
                           dimnames = list(
                             pars_to_sample$names,
                             pars_to_sample$names)),
  pars_obs = pars_obs, 
  pars_lprior = pars_lprior,
  model_params = model_params,
  sircovid_model = sircovid_model
)

expect_false(all(X2$results == X3$results))

})

test_that("pmcmc with new model", {
  data <- readRDS("hospital_model_data.rds")
  sircovid_model <- hospital_model()
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
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  pars_to_sample <- data.frame(
    names=c('beta_start', 'start_date'),
    init=c(0.14, sircovid_date("2020-02-07")),
    min=c(0, 0),
    max=c(1, 1e6),
    discrete=c(FALSE, TRUE),
    stringsAsFactors = FALSE)
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'start_date' = function(pars) 0)
  cmp <- readRDS("reference_pmcmc_hosp.rds")

  n_mcmc <- 10
  set.seed(1)

  X <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    pars_to_sample = pars_to_sample,
    pars_lprior = pars_lprior,
    proposal_kernel = matrix(c(0.001^2, 0,
                       0, 0.5^2),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(
                       c('beta_start', 'start_date'),
                       c('beta_start', 'start_date')))
  )

  expect_is(X, 'pmcmc')
  expect_setequal(names(X), c('inputs', 'results', 'states', 'acceptance_rate', 'ess'))
  expect_equal(dim(X$results), c(n_mcmc + 1L, 5))
  expect_equal(dim(X$states), c(n_mcmc + 1L, 289))
  expect_equivalent(X[-1], cmp[-1])
  
  
  # Test that the parameter order doesn't matter
  pars_to_sample <- data.frame(
    names=c('start_date', 'beta_start'),
    init=c(sircovid_date("2020-02-07"), 0.14),
    min=c(0, 0),
    max=c(1e6, 1),
    discrete=c(TRUE, FALSE),
    stringsAsFactors = FALSE)
  set.seed(1)
  X_reordered <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    pars_to_sample = pars_to_sample,
    pars_lprior = pars_lprior,
    proposal_kernel = matrix(c(0.001^2, 0,
                               0, 0.5^2),
                             nrow = 2, byrow = TRUE,
                             dimnames = list(
                               c('beta_start', 'start_date'),
                               c('beta_start', 'start_date')))
  )
  # Results will be different as random numbers generated in a different
  # order, but check that the right proposal has been made (0.5^2 vs 0.001^2 makes
  # a big difference)
  expect_true(all(X_reordered$ess[names(X_reordered$ess) != "log_prior"] > 0))
  expect_true(all(X_reordered$acceptance_rate[names(X_reordered$acceptance_rate) != "log_prior"] > 0))
  
  
  ## test generalised version
  
  pars_to_sample <- data.frame(
    names=c('beta_start',
            'beta_end', 
            'start_date',  
            'gamma_triage', 
            'gamma_hosp_R', 
            'gamma_hosp_D', 
            'gamma_ICU_R', 
            'gamma_ICU_D', 
            'gamma_stepdown'),
    init=c(0.14, 
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
          0),
    max=c(1,
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

  
  cmp <- readRDS("reference_pmcmc_gammas.rds")

  set.seed(2)
  X2 <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    pars_to_sample = pars_to_sample,
    pars_lprior = pars_lprior,
    proposal_kernel = proposal_kernel,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs
  )
  
  expect_is(X2, 'pmcmc')
  expect_setequal(names(X2), c('inputs', 'results', 'states', 'acceptance_rate', 'ess'))
  expect_equal(dim(X2$results), c(n_mcmc + 1L, length(pars_to_sample$names) + 3L))
  expect_equal(dim(X2$states), c(n_mcmc + 1L, 289))
  expect_equivalent(X2[-1], cmp[-1])
  
  # test that pars_obs can be modified
  pars_lprior <- list('beta_start'     = function(pars) log(1e-10),
                      'start_date'     = function(pars) 0,
                      'phi_general'    = function(pars) dnorm(pars['phi_general'], 0.95, 0.01))
  pars_to_sample <- data.frame(
    names=c('start_date', 'beta_start', 'phi_general'),
    init=c(sircovid_date("2020-02-07"), 0.14, 0.95),
    min=c(0, 0, 0),
    max=c(1e6, 1, 1),
    discrete=c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE)
  X_reordered <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    pars_to_sample = pars_to_sample,
    pars_lprior = pars_lprior,
    proposal_kernel = matrix(c(0.001^2, 0, 0,
                               0, 0.5^2, 0,
                               0, 0, 0.01^2),
                             nrow = 3, byrow = TRUE,
                             dimnames = list(
                               c('beta_start', 'start_date', 'phi_general'),
                               c('beta_start', 'start_date', 'phi_general')))
  )

})



test_that("pmcmc will run with multiple chains" , {
  
  data <- readRDS("hospital_model_data.rds")
  
  sircovid_model <- hospital_model()
  
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
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'beta_end' = function(pars) 0,
                     'start_date' = function(pars) 0)

  proposal_kernel <- matrix(c(0.001^2, 0, 0,
                             0, 0.001^2, 0,
                             0,       0, 0.5^2), 
                           nrow = 3, byrow = TRUE,
                           dimnames = list(
                             pars_to_sample$names,
                             pars_to_sample$names))

  n_mcmc <- 10
  n_chains <- 2
  set.seed(1)

  X <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    pars_to_sample = pars_to_sample,
    pars_lprior = pars_lprior,
    proposal_kernel = proposal_kernel,
    model_params = model_params,
    pars_obs = pars_obs, 
    n_chains = n_chains
  )

  expect_is(X, 'pmcmc_list')
  expect_equal(length(X$chains), n_chains)

  
  ## test create_master_chain
  Y <- create_master_chain(X, burn_in = 1L) 
  expect_error(create_master_chain(X, burn_in = 1e6), 
               'burn_in is greater than chain length')
  expect_error(create_master_chain(X, burn_in = '1'), 
               'burn_in must be an integer')
  expect_error(create_master_chain(X, burn_in = -1), 
               'burn_in must not be negative')
  expect_error(create_master_chain(X$chains, burn_in = 1L), 
               'x must be a pmcmc_list object')


  
  # Summary run, but not checked
  summary(X, burn_in = 1)
  ## plot called but not checked
  plot(X, burn_in = 1)
  
})

test_that("pmcmc error cases", {
  
  data <- read.csv(sircovid_file("extdata/example.csv"),
                   stringsAsFactors = FALSE)
  sircovid_model <- basic_model()
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'beta_end' = function(pars) 0,
                     'start_date' = function(pars) 0)
  proposal_kernel <- matrix(c(0.001^2, 0, 0,
                              0, 0.001^2, 0,
                              0,       0, 0.5^2), 
                            nrow = 3, byrow = TRUE,
                            dimnames = list(
                              pars_to_sample$names,
                              pars_to_sample$names))
  
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

  ## beta_start too low
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(-1, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    'initial parameters are outside of specified range'
  )
  
  ## start_date too late
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date(data$date[3])),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "start date must not be before first date of supplied data"
  )
  
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end'),
    init=c(0.14, 0.14*0.238),
    min=c(0, 0),
    max=c(1, 1),
    discrete=c(FALSE, FALSE),
    stringsAsFactors = FALSE) 
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "Turning off beta and start date sampling unsupported"
  )
  
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      pars_obs = pars_obs, 
      n_chains = n_chains,
      model_params = generate_parameters(
        sircovid_model = sircovid_model,
        transmission_model = "POLYMOD",
        beta = c(0.1, 0.1, 0.1),
        beta_times = sircovid_date(c("2020-02-02", "2020-03-01", "2020-04-01")),
        hosp_transmission = 0,
        ICU_transmission = 0,
        trans_profile = 1,
        trans_increase = 1,
        dt = 0.25)
    ),
    "Set beta variation through generate_beta_func in sircovid_model, not model_params"
  )

  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(-0.01, 0.03, sircovid_date("2020-02-01")),
    min=c(-1, -1, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE) 
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "beta_start must not be negative"
  )
  
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.123, -0.03, sircovid_date("2020-02-01")),
    min=c(-1, -1, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "beta_end must not be negative"
  )
  

 pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, as.character(sircovid_date(data$date[1]))),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "pars_min entries must be numeric"
  )
  
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, as.character(sircovid_date(data$date[1]))),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "pars_max entries must be numeric"
  )

  # incorrect format for proposal_kernel
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains,
      proposal_kernel = list('beta_start' = 0.5,
                     'beta_end' = 0.5,
                    'start_date' = 3),
    ),
    "proposal_kernel must be a matrix or vector with names corresponding to the parameters being sampled"
  )
  
  
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, 0.5, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "pars_discrete entries must be logical"
  )

  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains,
      output_proposals = 0:1
    ),
    "output_proposals must be either TRUE or FALSE"
  )
  
  
  
  ### checks on supplied log prior function

  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date'),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07")),
    min=c(0, 0, 0),
    max=c(1, 1, 1e6),
    discrete=c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains,
      pars_lprior = list('beta_start' = function(pars) {
        dunif(pars, min = 0, max = 1e6, log = TRUE)
      },
      'beta_end' = function(pars) 0,
      'start_date' = function(pars) 0)
    ),
    'log_prior must return a single numeric representing the log prior'
  )
  

  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains,
      pars_lprior = list('beta_start' = function(pars) {
        sum(dunif(pars, min = 1e6-1, max = 1e6, log = TRUE))
      },
      'beta_end' = function(pars) 0,
      'start_date' = function(pars) 0)
    ),
    'initial parameters are not compatible with supplied prior'
  )
  
  # Parameters not in pars_obs or model_params
  pars_to_sample <- data.frame(
    names=c('beta_start', 'beta_end', 'start_date', "amma_triage"),
    init=c(0.14, 0.14*0.238, sircovid_date("2020-02-07"), 0.5099579),
    min=c(0, 0, 0, 0),
    max=c(1, 1, 1e6, 1),
    discrete=c(FALSE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE)
  proposal_kernel <- matrix(c(0.001^2, 0, 0, 0,
                              0, 0.001^2, 0, 0,
                              0,       0, 0.5^2, 0,
                              0, 0, 0, 0.1^2), 
                            nrow = 4, byrow = TRUE,
                            dimnames = list(
                              pars_to_sample$names,
                              pars_to_sample$names))
  pars_lprior = list('beta_start' = function(pars) log(1e-10),
                     'beta_end' = function(pars) 0,
                     'start_date' = function(pars) 0,
                     'amma_triage' = function(pars) 0)
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_lprior = pars_lprior,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      n_chains = n_chains
    ),
    "Don't know how to update parameter: amma_triage"
  )
  
})

test_that("reflect_proposal", {
  
  expect_equal(object = reflect_proposal(x = 6, floor = 1, cap = 5), 
              expected = 4)
  
  expect_equal(object = reflect_proposal(x = 0, floor = 1, cap = 5), 
               expected = 2)
  
  expect_equal(object = reflect_proposal(x = 10, floor = 1, cap = 5), 
               expected = 2)
  

  
  # check that the function behaves as expected when passed a vector
  n <- 10
  tmp <- data.frame(x = rnorm(n, 1), 
                    floor = runif(n))
  tmp$cap <- tmp$floor + runif(n)
  
  X <- with(tmp, reflect_proposal(x, floor, cap))
  Y <- with(tmp, mapply(FUN = reflect_proposal, 
                        x = x, 
                        floor = floor, 
                        cap = cap))
  
  expect_equal(X, Y)
  
})

context("pmcmc")



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
    beta_times = '2020-01-01',
    hosp_transmission = 0,
    ICU_transmission = 0,
    trans_profile = 1,
    trans_increase = 1,
    dt = 1/4
  )

  X <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = list(phi_ICU = 0.95,
                    k_ICU = 2,
                    phi_death = 926 / 1019,
                    k_death = 2,
                    exp_noise = 1e6),
    pars_to_sample = c(
      'beta_start',
      'start_date'
    ),
    pars_init = list('beta_start' = 0.14,
                     'start_date' = as.Date("2020-02-07")),
    pars_min = list('beta_start' = 0,
                    'start_date' = 0),
    pars_max = list('beta_start' = 1,
                    'start_date' = 1e6),
    proposal_kernel = matrix(c(0.001^2, 0,
                       0, 0.5^2),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(
                       c('beta_start', 'start_date'),
                       c('beta_start', 'start_date'))),
    pars_discrete = list('beta_start' = FALSE,
                         'start_date' = TRUE)
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

  pars_to_sample <- c('beta_start', 'beta_end', 'start_date')
  pars_init <- list('beta_start' = 0.14, 
                   'beta_end' = 0.14*0.238,
                   'start_date' = as.Date("2020-02-07"))
  pars_min <- list('beta_start' = 0, 
                  'beta_end' = 0,
                  'start_date' = 0)
  pars_max <- list('beta_start' = 1, 
                  'beta_end' = 1,
                  'start_date' = 1e6)
  pars_discrete <- list('beta_start' = FALSE,
                       'beta_end' = FALSE,
                       'start_date' = TRUE)
  
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  
  
  set.seed(2)
 X2 <- pmcmc(
   data = data,
   n_mcmc = n_mcmc,
   pars_to_sample = pars_to_sample,
   pars_init = pars_init,
   pars_min = pars_min,
   pars_max = pars_max, 
   pars_discrete = pars_discrete,
   proposal_kernel = matrix(c(0.001^2, 0, 0,
                              0, 0.001^2, 0,
                              0,       0, 0.5^2), 
                            nrow = 3, byrow = TRUE,
                            dimnames = list(
                              pars_to_sample,
                              pars_to_sample)),
   pars_obs = pars_obs, 
   model_params = model_params,
   sircovid_model = sircovid_model
 )

 expect_equal(dim(X2$results), c(n_mcmc + 1L, 6))
 expect_equivalent(X2[-1], cmp[-1])

 
 ## set likelihood to accept every time with outlandish proposals
 set.seed(1)
 Y <- pmcmc(
   data = data,
   n_mcmc = n_mcmc,
   sircovid_model = sircovid_model,
   pars_to_sample = pars_to_sample,
   pars_init = pars_init,
   pars_min = pars_min,
   pars_max = pars_max, 
   pars_discrete = pars_discrete,
   model_params = model_params,
   pars_obs = pars_obs,
   proposal_kernel = matrix(c(1e2, 0, 0,
                      0, 1e2, 0,
                      0,  0, 1e4),
                    nrow = 3, byrow = TRUE,
                    dimnames = list(
                      c('beta_start', 'beta_end', 'start_date'),
                      c('beta_start', 'beta_end', 'start_date')
                    )),
   log_likelihood = function(pars, ...) {
     list('log_likelihood' = 0,
          'sample_state' = rep(1, 238))
   }
 )
 expect_equal(Y$results$log_likelihood, rep(0, n_mcmc + 1L))
 # check that all proposals have been accepted
 expect_true(all(diff(Y$results$beta_start) != 0))
 # check that all start_dates are before data
 expect_true(max(Y$results$start_date) <= data$date[1])
 # check beta_start is in specified range
 expect_true(min(Y$results$beta_start) > 0)
 expect_true(max(Y$results$beta_start) < 1)
 expect_true(min(Y$results$beta_end) > 0)
 expect_true(max(Y$results$beta_end) < 1)

 ## check that proposing jumps of size zero results in the initial parameter being retained

 Z <- pmcmc(
   data = data,
   n_mcmc = n_mcmc,
   pars_to_sample = pars_to_sample,
   pars_init = pars_init,
   pars_min = pars_min,
   pars_max = pars_max, 
   pars_discrete = pars_discrete,
   sircovid_model = sircovid_model,
   pars_obs = pars_obs,
   model_params = model_params,
   proposal_kernel = matrix(rep(0, 9),
                    nrow = 3, byrow = TRUE,
                    dimnames = list( pars_to_sample, pars_to_sample)),
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
  pars_init = pars_init,
  pars_min = pars_min,
  pars_max = pars_max, 
  pars_discrete = pars_discrete,
  sircovid_model = sircovid_model,
  model_params = model_params,
  pars_obs = list(phi_ICU = 0.95,
                  k_ICU = 2,
                  phi_death = 926 / 1019,
                  k_death = 2,
                  exp_noise = 1e6),
  proposal_kernel = matrix(c( 0.001^2,       0, 0.001*0.5*0.6,
                                    0, 0.001^2, 0,
                        0.001*0.5*0.6,       0, 0.5^2),
                   nrow = 3, byrow = TRUE,
                   dimnames = list( pars_to_sample,pars_to_sample))
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
    beta_times = '2020-01-01',
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
  cmp <- readRDS("reference_pmcmc_hosp.rds")

  n_mcmc <- 10
  set.seed(1)

  X <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    model_params = model_params,
    pars_obs = pars_obs,
    pars_to_sample = c(
      'beta_start',
      'start_date'
    ),
    pars_init = list('beta_start' = 0.14,
                     'start_date' = as.Date("2020-02-07")),
    pars_min = list('beta_start' = 0,
                    'start_date' = 0),
    pars_max = list('beta_start' = 1,
                    'start_date' = 1e6),
    proposal_kernel = matrix(c(0.001^2, 0,
                       0, 0.5^2),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(
                       c('beta_start', 'start_date'),
                       c('beta_start', 'start_date'))),
    pars_discrete = list('beta_start' = FALSE,
                         'start_date' = TRUE)
  )

  expect_is(X, 'pmcmc')
  expect_setequal(names(X), c('inputs', 'results', 'states', 'acceptance_rate', 'ess'))
  expect_equal(dim(X$results), c(n_mcmc + 1L, 5))
  expect_equal(dim(X$states), c(n_mcmc + 1L, 289))
  expect_equivalent(X[-1], cmp[-1])
  
  
  ## test generalised version
  
  pars_to_sample <- c('beta_start','beta_end', 'start_date',  'gamma_triage', 'gamma_hosp_R', 
                      'gamma_hosp_D', 'gamma_ICU_R', 'gamma_ICU_D', 'gamma_stepdown')
  
  proposal_kernel <- diag(length(pars_to_sample)) * 0.01^2
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- pars_to_sample
  proposal_kernel['start_date', 'start_date'] <- 25
  
  
  cmp <- readRDS("reference_pmcmc_gammas.rds")

  set.seed(2)
  X2 <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    pars_to_sample = pars_to_sample,
    proposal_kernel = proposal_kernel,
    sircovid_model = sircovid_model
  )
  
  expect_is(X2, 'pmcmc')
  expect_setequal(names(X2), c('inputs', 'results', 'states', 'acceptance_rate', 'ess'))
  expect_equal(dim(X2$results), c(n_mcmc + 1L, length(pars_to_sample) + 3L))
  expect_equal(dim(X2$states), c(n_mcmc + 1L, 289))
  expect_equivalent(X2[-1], cmp[-1])
  
  

})

test_that("pmcmc will run with multiple chains" , {
  
  data <- readRDS("hospital_model_data.rds")
  
  sircovid_model <- hospital_model()
  
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
  pars_obs <- list(
    phi_general = 0.95,
    k_general = 2,
    phi_ICU = 0.95,
    k_ICU = 2,
    phi_death = 926 / 1019,
    k_death = 2,
    exp_noise = 1e6
  )
  cmp <- readRDS("reference_pmcmc_hosp.rds")
  
  pars_to_sample <- c('beta_start', 'beta_end', 'start_date')
  pars_init <- list('beta_start' = 0.14, 
                    'beta_end' = 0.14*0.238,
                    'start_date' = as.Date("2020-02-07"))
  pars_min <- list('beta_start' = 0, 
                   'beta_end' = 0,
                   'start_date' = 0)
  pars_max <- list('beta_start' = 1, 
                   'beta_end' = 1,
                   'start_date' = 1e6)
  pars_discrete <- list('beta_start' = FALSE,
                        'beta_end' = FALSE,
                        'start_date' = TRUE)
  proposal_kernel <- matrix(c(0.001^2, 0, 0,
                             0, 0.001^2, 0,
                             0,       0, 0.5^2), 
                           nrow = 3, byrow = TRUE,
                           dimnames = list(
                             pars_to_sample,
                             pars_to_sample))

  n_mcmc <- 10
  n_chains <- 2
  set.seed(1)

  X <- pmcmc(
    data = data,
    n_mcmc = n_mcmc,
    sircovid_model = sircovid_model,
    pars_to_sample = pars_to_sample,
    pars_init = pars_init,
    pars_min = pars_min,
    pars_max = pars_max, 
    pars_discrete = pars_discrete,
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
  
  pars_to_sample <- c('beta_start', 'beta_end', 'start_date')
  pars_init <- list('beta_start' = 0.14, 
                    'beta_end' = 0.14*0.238,
                    'start_date' = as.Date("2020-02-07"))
  pars_min <- list('beta_start' = 0, 
                   'beta_end' = 0,
                   'start_date' = 0)
  pars_max <- list('beta_start' = 1, 
                   'beta_end' = 1,
                   'start_date' = 1e6)
  pars_discrete <- list('beta_start' = FALSE,
                        'beta_end' = FALSE,
                        'start_date' = TRUE)
  proposal_kernel <- matrix(c(0.001^2, 0, 0,
                              0, 0.001^2, 0,
                              0,       0, 0.5^2), 
                            nrow = 3, byrow = TRUE,
                            dimnames = list(
                              pars_to_sample,
                              pars_to_sample))
  
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

  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_init = list(0.5,
                       as.Date("2020-02-01")),
      pars_to_sample = pars_to_sample,
      pars_min = pars_min,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
      
    ),
    "pars_init must be a named list corresponding to the parameters being sampled"
  )

  ## beta_start too low
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_obs = pars_obs, 
      model_params = model_params,
      pars_init = list('beta_start' = -1,
                       'beta_end' = 0,
                       'start_date' = as.Date("2020-02-01")),
      pars_to_sample = pars_to_sample,
      pars_min = pars_min,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    'initial parameters are outside of specified range'
  )
  
  ## start_date too late
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_obs = pars_obs, 
      model_params = model_params,
      pars_init = list('beta_start' = 0.123,
                       'beta_end' = 0.03,
                       'start_date' = as.Date(data$date[3])),
      pars_to_sample = pars_to_sample,
      pars_min = pars_min,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    'initial parameters are outside of specified range'
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_obs = pars_obs, 
      model_params = model_params,
      pars_init = list('beta_start' = 0.123,
                       'beta_end' = 0.03,
                       'start_date' = as.Date(data$date[1])), 
      pars_min = list('beta_start' = 0, 
                      'beta_end' = 0,
                      'start_date' = 0),
      pars_to_sample = pars_to_sample,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "'start_date' must be less than the first date in data"
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_obs = pars_obs, 
      model_params = model_params,
      pars_to_sample = c(
        'beta_start',
        'beta_end'
      ),
      pars_init = pars_init, 
      pars_min = pars_min,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "Turning off beta and start date sampling unsupported"
  )
  
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_to_sample = pars_to_sample,
      pars_init = pars_init,
      pars_min = pars_min,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel,
      pars_obs = pars_obs, 
      model_params = generate_parameters(
        sircovid_model = sircovid_model,
        transmission_model = "POLYMOD",
        beta = c(0.1, 0.1, 0.1),
        beta_times = c("2020-02-02", "2020-03-01", "2020-04-01"),
        hosp_transmission = 0,
        ICU_transmission = 0,
        trans_profile = 1,
        trans_increase = 1,
        dt = 0.25)
    ),
    "Set beta variation through generate_beta_func in sircovid_model, not model_params"
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_init = list('beta_start' = -0.01,
                       'beta_end' = 0.03,
                       'start_date' =  as.Date("2020-02-01")), 
      pars_min = list('beta_start' = -1 ,
                      'beta_end' = -1,
                      'start_date' = 0),
      pars_to_sample = pars_to_sample,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "beta_start must not be negative"
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_init = list('beta_start' = 0.123,
                       'beta_end' = -0.03,
                       'start_date' =  as.Date("2020-02-01")), 
      pars_min = list('beta_start' = -1 ,
                      'beta_end' = -1,
                      'start_date' = 0),
      pars_to_sample = pars_to_sample,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "beta_end must not be negative"
  )
  
  
  # incorrect names supplied to pars_min
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_min  = list(0.5, 0.03,0),
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs
    ),
    "pars_min must be a named list corresponding to the parameters being sampled"
  )
  
  # incorrect format for pars_min
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_min  = c('beta_start' = 0.5,
                    'beta_end' = 0.03, 
                       'start_date' = 0), 
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "pars_min must be a named list corresponding to the parameters being sampled"
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_min  = list('beta_start' = 0.5,
                       'beta_end' = 0.03, 
                    'start_date' = as.Date(data$date[1])),
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_max = pars_max, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "pars_min entries must be numeric"
  )
  
  # incorrect names supplied to pars_max
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_max  = list(0.5, 1, 1e3), 
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs
    ),
    "pars_max must be a named list corresponding to the parameters being sampled"
  )
  
  # incorrect format for pars_max
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_max  = c('beta_start' = 0.5,
                    'beta_end' = 0.03, 
                    'start_date' = 1e3),
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "pars_max must be a named list corresponding to the parameters being sampled"
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_max  = list('beta_start' = 0.5,
                       'beta_end' = 1, 
                       'start_date' = as.Date(data$date[1])),
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_discrete = pars_discrete,
      proposal_kernel = proposal_kernel
    ),
    "pars_max entries must be numeric"
  )
  
  
  # incorrect names supplied to proposal_kernel
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      proposal_kernel = list(0.5, 0.5, 3),
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs
    ),
    "proposal_kernel must be a matrix or vector with names corresponding to the parameters being sampled"
  )
  
  # incorrect format for proposal_kernel
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      model_params = model_params,
      pars_obs = pars_obs, 
      proposal_kernel = list('beta_start' = 0.5,
                     'beta_end' = 0.5,
                    'start_date' = 3),
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      pars_discrete = pars_discrete
    ),
    "proposal_kernel must be a matrix or vector with names corresponding to the parameters being sampled"
  )
  
  
  # incorrect names supplied to pars_discrete
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_discrete = list(FALSE, FALSE, TRUE)
    ),
    "pars_discrete must be a named list corresponding to the parameters being sampled"
  )
  # incorrect format for pars_discrete
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_discrete = c('beta_start' = FALSE, 
                        'beta_end' = FALSE,
                        'start_date' = TRUE)
    ),
    "pars_discrete must be a named list corresponding to the parameters being sampled"
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      model_params = model_params,
      pars_obs = pars_obs, 
      pars_discrete  = list('beta_start' = 0.5,
                            'beta_end' = 0.5,
                       'start_date' = as.Date(data$date[1]))
    ),
    "pars_discrete entries must be logical"
  )

  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      output_proposals = 0:1
    ),
    "output_proposals must be either TRUE or FALSE"
  )
  
  
  
  ### checks on supplied log prior function

  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_prior = function(pars) {
        dunif(pars, min = 0, max = 1e6, log = TRUE)
      }
    ),
    'log_prior must return a single numeric representing the log prior'
  )
  

  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_prior = function(pars) {
        sum(dunif(pars, min = 1e6-1, max = 1e6, log = TRUE))
      }
    ),
    'initial parameters are not compatible with supplied prior'
  )
  
  # checks on supplied log likelihood function
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars) {
        dunif(pars, min = 0, max = 1e6, log = TRUE)
      }
    ),
    'log_likelihood function must be able to take unnamed arguments'
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars, ...) {
        x <- sum(dunif(pars, min = 0, max = 1e6, log = TRUE))
        list('log_likelihood' = x)
      }
    ),
    'log_likelihood function must return a list containing elements log_likelihood and sample_state'
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars, ...) {
        x <- sum(dunif(pars, min = 0, max = 1e6, log = TRUE))
        list('log_likelihood' = x, "sample_state" = x)
      }
    ),
    'sample_state must be a vector of non-negative numbers'
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars, ...) {
        x <- sum(dunif(pars, min = 0, max = 1e6, log = TRUE))
        list(x, rep(1, 236))
      }
    ),
    'log_likelihood function must return a list containing elements log_likelihood and sample_state'
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars, ...) {
        x <- dunif(pars, min = 0, max = 1e6, log = TRUE)
        list('log_likelihood'= x, 'sample_state' = rep(1, 236))
      }
    ),
    'log_likelihood must be a single numeric representing the estimated log likelihood'
  )
  
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars, ...) {
        x <- sum(dunif(pars, min = 0, max = 1e6, log = FALSE))
        list('log_likelihood'= x, 'sample_state' = rep(1, 236))
      }
    ),
    'log_likelihood must be negative or zero'
  )
  
  expect_error(
    pmcmc(
      data = data,
      n_mcmc = n_mcmc,
      sircovid_model = sircovid_model,
      pars_init = pars_init,
      pars_to_sample = pars_to_sample,
      pars_min = pars_min, 
      pars_max = pars_max,
      proposal_kernel = proposal_kernel,
      pars_discrete = pars_discrete,
      model_params = model_params,
      pars_obs = pars_obs, 
      log_likelihood = function(pars, ...) {
        x <- sum(dunif(pars, min = 0, max = 1e6, log = TRUE))
        list('log_likelihood'= x, 'sample_state' = rep(-1, 236))
      }
    ),
    'sample_state must be a vector of non-negative numbers'
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
